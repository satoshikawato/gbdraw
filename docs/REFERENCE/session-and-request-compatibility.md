[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md)

# Session and request compatibility

Current writers emit session version 49 and canonical `renderRequest` schema 9.
Session 49 lists the drawings of a project in `drawings`. The Web app and Python
read a Session 46–48 into drawings themselves; a Session 27–44 is converted by
Python, once, at Load (the Web app waits for it); a Session 45 is not read.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 49 | 27–33, 39–42, 44, and 46–49 |
| Canonical `renderRequest` | 9 | 1, 2, 5, 6, 7, 8, and 9 |
| Web file bindings | 2 | 1; 2 in sessions 41–42, 44, and 46–48 |

Sessions 48 and 49 store the saved protein results (LOSATP tables, Similarity groups,
Collinear blocks) of a reverse-complemented Linear record in the record's
search frame, as it stores comparison tables. Sessions 46 and 47 stored them
after the reverse complement. A replay of an older Session draws as before,
and a re-save or `upgrade_session_document()` writes Session 49. A saved result
that matches neither strand of its record is still drawn, but the Session
cannot be saved again (`SessionConversionError`); generate the comparison
again.

Session 48 also records which rows, colors, and Blacklist words each file card
of the Web app loaded, so **Remove** on the card takes back only that file's
effect. In an older Session, Whitelist and Priority rows read as rows added by
hand, and Remove keeps them; Specific rules rows that came from a rules file
were marked so already, and Remove takes them back. Rows that a built-in
**Preset scheme** added count as yours, and Remove keeps them; an older Session
stored them like a rules file's rows, so there Remove takes them back. In a Session 48 that the
CLI or the Python API writes, the settings read from its color, priority,
and whitelist files count as those files' cards. `--label_blacklist` takes
keywords, so a CLI Session records it as keywords only; a CLI Session written
by gbdraw 0.13.0 stored an existing path given to it as a Blacklist File, and a
CLI replay keeps that card.

A Result saved by gbdraw 0.12.0 to 0.13.0, or by the web app from 2026-06-18 to
2026-07-19 (Sessions 27 to 32), names each
comparison match by `data-gbdraw-pairwise-match-id`. The web app renames it to
`data-gbdraw-match-id`, the only name gbdraw writes now, when it loads such a
Session, so match popups and selection work as before.

The palette `herbal_garden` was named `hearbal_garden` through gbdraw 0.14.2.
`-p/--palette`, `FeatureOptions(palette=...)`, `ColorOptions`, and a Session or
request that names `hearbal_garden` draw `herbal_garden`, and a Session that the
CLI or the Python API saves names `herbal_garden`; a CLI Session keeps the
command line as typed. The web app loads the drawing's palette and a queued
palette as `herbal_garden`. Until the next **Generate**, a Session it saves
keeps the former name in the palette of the shown Result and in its render
request; both still draw `herbal_garden`.

Session 44 documents written with canonical request schema 7 or 8 remain
readable. Current saves promote that request to schema 9 and the Session to 48,
and keep the saved preview until the next Generate. Schema 9 adds
`diagramOptions.featureOverrides` (see
[Feature identity overrides](typed-requests.md#feature-identity-overrides));
a schema-8 request reads as one without overrides.

The Python reader (`gbdraw ... --session` and the Python API) checks a
schema-9 request against `describe_canonical_request()`: it refuses a key
outside it, such as an older writer's `depthTrack*` input, and a value of
another kind, such as `true` or `false` where a number belongs. At every
request schema it refuses `conservationTableFile` and `comparisonTableFile`
(a table is stored as a resource); anything but text (or `null` where a
field is optional) at an annotation id, label, legend label, selector key or
value, fill,
`generatedProteinComparison.mode`, `conservationReference` or
`pairwiseMatchStyle`; and anything but `true` or `false` at `wrapsOrigin` or
hatch `cross`. From schema 9 it also refuses another case, outer spaces or
`""` in these three names, and the five comparison mode settings
(`collinearityUnitMode`, `collinearityAnchorMode`, `collinearitySearchScope`,
`collinearityColorMode`, `orthogroupMembershipMode`) hold one of their names;
older schemas also read another case, outer spaces, `-` for `_`, an alias, or
`""`. An annotation `lane`, `start` or `end` may hold a fraction, which reads
as its integer part. Writers store a track slot's `renderer` and `side` as
names, and the reader still reads another spelling. For a form listed here
that earlier versions read, the reader's message names the refused path. The
Web app still loads a Session that holds such a form; until the Web repairs it
at Load, Generate and any re-render of the stored request (label reflow,
alignment, record transform) show "An input value is invalid"
(`INPUT_INVALID`, field `schema`).

Session 49 holds a project and its list of drawings. The project keeps
`format`, `version`, `createdAt`, `title`, `resources`, `ui` (with
`ui.activeDrawingId`, the drawing that opens), and `legacyArtifacts`. Each entry
of `drawings` is one drawing in tab order: its `id`, `name`, and `mode`
(`circular` or `linear`), and its own request, Results, settings (`config`),
per-feature edits (`features`), editor state, view (`ui`), file bindings
(`webFiles`), and saved LOSAT results. A project may hold several drawings of one
mode. The top-level `modes`, `otherModeResult`, `renderRequest`, `results`,
`config`, and `features` of Sessions 46 to 48 are not part of Session 49.
Reading a Session 46 to 48 moves each mode's part into its own drawing.

Sessions 46 to 48 keep the Web app's settings and edits of each diagram mode in its
own slice, `modes.circular` and `modes.linear`:

- `config`: the draft configuration.
- `features`: the per-feature edits.
- `editorState`: the Legend edits and per-feature strokes.
- `ui`: the mode's layout, padding, and pending palette.

A field missing from a slice reads as that mode's default. A Session 46 to 48 rejects
the top-level `config` and `features` of Sessions 44 and older. Reading an
older Session splits its draft into the two slices. A Session that the CLI or
the Python API writes from a request has no slices; a CLI re-save keeps or
splits the source's draft.

Sessions 47 and 48 store each slice's Legend color and stroke edits
(`editorState.legend.colorOverrides` and `strokeOverrides`) under the row's
key: the caption gbdraw drew for the row (`originalCaption`), not the name the
row shows. Two rows can then show one name and keep their own edits. Session 46
stored these edits under the shown name. Reading a Session 46, or an older
Session, gives each row the edit that Session 46 applied to it. An edit stored
under a renamed row's shown name moves to that row's key. When several rows
show that name, the first listed row takes the edit, and a deleted row only
when no listed row shows it. An edit stored under the original caption of a
renamed row that no row shows any more stays on that row, where Session 46
applied it. When two edits reach one row, the one Session 46 applied last is
kept; a stroke takes its color and its width each from the last edit that sets
it. The Web app, a CLI re-save, and `upgrade_session_document()` read
these edits the same way.

Sessions 46 to 48 store the per-feature edits (Feature visibility, Label visibility,
and label text) in each slice's `features.featureOverrides`. There is one row
per original-source feature, named by `recordKey` and `biologicalFeatureId`, as
in the request. Feature placement drafts and the `featureIdentity` targets of
annotations made from selected features are kept per mode the same way. Both
modes can use the same record key, so a request carries only the rows and
targets of its own mode. A row keeps the record key of the view it was last
drawn in; the next Generate that draws the record under another key moves the
row to it, so a Session saved by 0.14.2 or earlier keeps its edits when the
view changes. The feature catalog is schema 5. It records the
selector values each drawn feature had (`drawnSelector`).

When you load an older Session, the Web app moves its per-feature edits onto
the features they name:

- A Session 44, or a Web Session 31–33, keys its edits by rendered feature ID.
  A rendered ID in the saved feature catalog names its feature. Otherwise a
  copy or record suffix is removed when exactly one feature remains.
- A Session before 40 is matched through its GenBank sources, read again with
  their crops and orientations (the drawn hash and record position of each
  rendered ID). Without readable sources, it is matched through its saved
  feature metadata, for records drawn without a crop or reverse complement.
- If the Web app's diagram runtime cannot start or fails while reading the
  sources, the Load fails and keeps the current Session.
- A CLI re-save and `upgrade_session_document()` read the sources again the
  same way. When a source cannot be read, they use the saved feature metadata,
  without the Web app's fallback to the features of the saved SVG.
- An edit that names no feature is dropped, and the Web app lists how many
  before it loads the Session.
- An older Feature visibility edit that hid every feature with its hash (each
  copy of a duplicated record) now applies to the edited feature only; the Web
  app lists how many before it loads the Session.
- Each moved edit goes to the slice of the mode of the Session's diagram. A
  Feature placement draft of a Session 41–44 goes to the slice of its lane
  side's mode. A Main placement, which reached requests of both modes, goes to
  both slices.
- An annotation target `hash=<hash>` of a Session 40–44 moves to a
  `featureIdentity` target in the Session's mode only when its record is drawn
  without a crop, reverse complement, or rotation and the hash names exactly
  one feature of the saved catalog, in that record. The Web app reports how
  many moved, and every other target loads unchanged.
- Since 0.14.0, `hash=` names a feature by its hash in its source record.
  Sessions up to 44 matched it against the drawn feature. On a cropped or
  reverse-complemented record, their Web app wrote drawn hashes in two places:
  the `hash` specific color rules of **This feature only** (release 0.13.0 wrote
  the rendered ID, with its `_record_<n>` suffix), which the request's color
  table repeats, and the `hash` rows of the request's label override table,
  which name a label text edit of a feature without a unique qualifier by its
  rendered ID.
  - Load and `upgrade_session_document()` give each such color rule the source
    hash of the features the saved diagram drew with it, through the saved
    catalog or the sources read again. `--session` replay does the same for the
    color and label override tables of the request.
  - A value whose features have its own hash as source hash, on a record drawn
    without a crop or reverse complement, stays as saved.
  - A value that names features of several source hashes is kept; the Web
    app lists how many before it loads, and the CLI reports how many. So is a color rule of a record drawn
    cropped or reverse-complemented that the CLI cannot read: the Session that
    `--session_output` writes for a `--session` replay does not read the
    sources of a Session 31–39 again, and no CLI path reads those of a Session
    27–30.
  - The Feature visibility rows, the draft's label override rows, and the
    `hash=` annotation targets of these Sessions hold source hashes or values
    a user typed. They load as saved and now match the feature with that
    source hash on every record.
  - A label whitelist row that the former **Whitelist only** action wrote for a feature without
    a unique qualifier holds its rendered ID and is not renamed either.
- A catalog of schema 3 or 4 reads as schema 5 without selector values until
  the next Generate. A feature whose rendered ID carries its source hash was
  drawn with its source coordinates, so its source values serve until then.

A CLI Circular batch or grid Session stores a cropped or reverse-complemented
record as drawn, so its `--session` replay gives those features the hashes of
the drawn copy, not of the original source record.

In a Session of any version, the Web app reads a Legend entry color with the
[Default colors](input-formats-and-tsv-schemas.md#styling-tables) forms. An
entry whose color is outside them is dropped, and the Load dialog names it.

Legend names are numbered when the legend is drawn, by the rule under
[Feature presentation](palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation),
and two readings keep older files usable:

- A specific color rule's caption that ends in its own rule's color ID, such
  as `Packaging [#ffe700]` on a rule colored `#ffe700`, is read without that
  ID, from a saved Session and from a color table file, also when you typed
  that name yourself. The Web app, the CLI `--session` replay, and the Python
  Session readers read it so in a Session of any version.
- Release 0.13.0 drew two tracks with one legend label as `Coverage` and
  `Coverage (2)`. Its Sessions load, and their legend row edits stay on the
  rows of the same tracks. The rows are now drawn `Coverage (1)` and
  `Coverage (2)` when their colors differ. With equal colors they are one
  `Coverage` row, and the edits saved for the second row are not used. Two
  GC content tracks always have equal colors.

A Session 46 to 48 can also hold `otherModeResult`: the Result set of the diagram
mode that the top-level set does not draw (Circular or Linear), written when
both modes have a Result. Its `renderRequest`, `results`, `editorState`
(`featureCatalog`, `alignmentResetReceipt`, and the Result's generated Legend
order, colors, and stroke defaults), `ui`, `runMetadata`, and `cliInvocation`
mirror the top-level fields; its request names resources in the top-level
`resources` table. It needs a top-level request of the other mode and at least
one Result, so a settings-only Session cannot hold it.

## Loading a file that cannot be fully restored

The Web app reads a Session or an older configuration file as far as it can.
Before any setting, Result, input file, History entry or screen changes, it
shows one **Load this file?** dialog that lists the values it will drop, replace
with a default, not use, or read with another meaning:

- settings their reader cannot use, named by their place in the file and the
  value used, for example `modes.circular.config.adv.depth_tracks[0].height`.
  A label, file name, ID or other free text is named by its place only. A
  configuration file loads over the current drawing, so a settings group
  (`form`, `adv`, `losat` or `circularConservation`) that is not an object, a
  value of another kind inside one, and `rules` or `whitelist` that is not a
  list (such as `"form": 5`) are not used: the drawing keeps its value, and
  the row names that value or says that the drawing keeps its current value;
- Legend positions outside the diagram mode in the saved layout preferences
  or a saved Result, with the position used;
- Legend entries with a color the app does not accept or that repeat an
  earlier entry, and Legend edits that a later edit or a renamed row replaces;
- per-feature edits and `hash=` values as listed above;
- joined or dropped table rows;
- Depth rows of the drawing's track stack that name a Depth series the drawing
  does not have; they are turned off and no longer name it. A track stack
  that cannot be read is read as the default stack. A slot value that cannot
  be used is left unset or read as its default, and an axis index outside the
  stack is read as the nearest place in it;
- Legend and stroke edits, original Legend colors and the Legend order of a
  saved Result that cannot be read;
- input file entries (a binding or a stored file) that name a file the Session
  does not hold or cannot read. A comparison sequence file or a Depth file
  keeps its place empty, so no later file changes series or track. A Circular
  comparison series whose BLAST table cannot be read is left out with its
  sequence files, with a row that names the series; the other series keep
  their own files and names (a Session 27–30 draft without series names its
  rings by `labels`, by place, so the left-out ring's label goes with it). A
  Linear comparison whose BLAST TSV cannot be read loads without that file,
  as after removing the file from its pair card;
- parts of the Session file that Load cannot read, each named by its place:
  the saved Results, the feature list and edits of the saved diagram, the
  stored input files, the saved diagram of the other mode, or the settings of
  a drawing; and fields that a Session of that version never held, such as
  `modes` in a Session older than version 46 or an unknown top-level key;
- feature information or comparison match sequences that cannot be restored.

Rows of per-feature edits, Feature visibility edits, `hash=` values, feature
information and comparison match sequences give a count or a statement, not a
place in the file. Legend entry rows name the entry, and table rows name the
table and its lines. A Session of version 44 or older keeps a setting that
both diagram modes read in one place, so such a setting it cannot use is one
row (`config.adv.label_rendering`); a Session of version 46 or later keeps each
mode's settings in its own place, and each place is a row
(`modes.circular.config.adv.label_rendering`, `modes.linear.config.adv.label_rendering`).

The dialog shows at most 12 lines, the last one with the number of the rest,
and says that Save Session afterwards writes the current format. Load reads
each version of a Session against the fields and values that version can hold,
and lists every value it drops or replaces. One reading still keeps a value as
stored: an unknown text value of `track_type`, `scale_style`, `labels_mode` or
`show_labels_linear`, in a Session or a configuration file, because the
request format gives them no fixed list of values to check against.

- **Load** loads the file. **Cancel** changes nothing. The file itself is
  never changed.
- No dialog opens for a reading that keeps the meaning (an older format,
  capitals, outer spaces, an alias, a number stored as its text, or a blank
  number, such as a cleared number field, which reads as the field's default
  in a Session and in a configuration file), a name that an earlier gbdraw
  wrote and that the Session now keeps elsewhere, or a value that had no
  effect on the figure: a cache that the next comparison or Generate computes
  again (the LOSAT caches, the protein identity manifest), run notices and
  annotation warnings that no longer match the saved Results, an empty value
  (`null`, `[]` or `{}`) of an optional part, a stored input file that nothing
  names, a Session 46 Legend edit that a later edit of the
  same row overwrote, a pending palette that Instant Preview drops, a name or
  description edit of a similarity group the saved Result does not have, a
  setting of the other diagram mode, view state such as zoom, or a value
  another stored value replaces. A file that loads without a loss, such as
  every Gallery Session, opens none.
- A Depth series without a Depth file loads as saved, and so does a Depth row
  the app keeps when you remove a Depth file (a row with its own size or
  position). Save Session before the first Generate of a mode leaves such an
  enabled row out of the saved render request, as it draws nothing; the
  Session's settings keep the row, and Generate still stops on its row issue.
- The **Replace the current work?** question, asked when there is unsaved
  work, still comes first, before you choose the file. The question about a
  Session font that is not available (see [Fonts](#fonts)) comes after
  this dialog, so **Cancel** here loads no font either. The notice of `hash=` annotation targets moved to
  the feature they name, the notice of a loaded configuration file, and the
  **Comparisons** card notices still follow the load.

### Saved diagram settings that cannot be read

Load reads the saved diagram's settings (the stored render request of request
schemas 5–9 and the other diagram mode's saved result) through one reader and
checks them against the request description. What it cannot use is repaired
instead of refused, and each repair is a row in the dialog:

- an unknown setting is dropped;
- an invalid value is reset to its default or left unset;
- an unreadable entry of a list is dropped;
- a part that cannot be read is left out with the settings inside it, or
  replaced by its default when it has one;
- a reference to a file or record the Session does not hold is left out. A
  saved comparison that the **Comparisons** card keeps read-only or asks
  about keeps its record references, which that card reports.

A font size of 0 or less in a Session (feature labels, Legend, scale and
ruler labels, ticks, definition text, plot title, GC content and Depth ticks,
or a Linear Record Labels style such as Name / Species) reads as Auto, the
default size, with the row `Settings: the value "0" (<place>) cannot be used;
its default is used.` `gbdraw render --session`, `gbdraw circular --session`,
`gbdraw linear --session` and the Python API read it the same way and print
that row as a `WARNING:` line. **Generate**, the CLI size options and the
Python options still refuse such a size.

When the saved diagram's settings cannot be read at all (for example, their
`schema` is not an integer), or are missing, the file loads as settings only if it holds the web app's own settings: a Session
of version 42, 44, 46, 47 or 48 saved by the web app. A row says that the saved
Results, the GenBank, GFF and FASTA inputs, the feature list and the caches
are left out. A Session of version 42, 44, 47 or 48 written by the CLI or the
Python API, and a Session of version 31–41, has no saved settings that can be
read without them, so it is refused; so is a Session of version 31–41 whose
request names a record file the Session does not hold. A `null` request
reads as a missing one in a file that holds saved Results or a feature list;
in a file without them it marks a settings-only Session. When the settings of the
other diagram mode's saved result cannot be read, that diagram and its Results
are left out, with a row.

The five comparison settings (unit, anchor, search scope, color and
membership) read the same spellings in the web app and in Python in request
schemas 5–8. In a schema 9 request the web app also reads another spelling of
a value (such as `top1` or ` ALL `) as that value, without a row, while the CLI
and the Python API refuse it. Eight search-scope spellings that only the web
app accepted (`adjacent_pairs`, `adjacent_pair`, `neighbor_pairs`,
`neighboring_pairs`, `all_records`, `all_pairs`, `all_record_pairs`, `global`)
are no longer read: in the saved diagram's settings, a generated comparison
that holds one is dropped, with a row; in a Session's comparison settings, the
default is used, with a row. A file that holds `top1` or `legacy` now loads.

When a part of a Session other than its settings cannot be read, Load drops
that part and lists it. A Result that stays but whose feature information
cannot be restored keeps its SVG preview, and **Generate** restores the
feature information. Only the cases below stop a Load.

The Web app refuses these files and shows an error instead:

- A file that is not JSON, including a broken gzip file and text that is not
  UTF-8.
- A file that is neither a gbdraw Session nor a configuration file, including
  a Session without its `format` or `version` and request schemas 3 and 4.
- A Session of version 26 or older (gbdraw 0.11.0 or older, or gbdraw.app
  before 0.12.0). The error says that an older gbdraw saved it.
- A Session of a version this reader does not read, such as the
  development-only Sessions 34–38, 43 and 45. The error names the version and
  the versions this gbdraw reads.
- A Session version or request schema newer than this reader. For a Session,
  the error names its version and the newest version this gbdraw reads.
- A file larger than 200 MB, or larger than 512 MB after decompression.
- A gzip file opened in a browser without gzip decompression, or any file in
  a browser without module Workers. A browser feature or Worker that fails
  during the Load also stops it with an error.
- A file with an unsafe name (`__proto__`, `constructor`, `prototype`).
- A file from which Load restores nothing: every part of its content that it
  holds (settings, the saved diagram, Results, input files, edits or caches)
  would be dropped or replaced by its default. A settings-only Session whose
  `modes` is not an object, or whose settings hold nothing or were all
  dropped, is such a file, and so is a Session of version 31–41, or one of
  version 42, 44, 47 or 48 written by the CLI or the Python API, whose saved
  settings cannot be read. A Session of version 27–30 that holds only its
  `format` and `version` loads.

Load refuses these before it asks anything or changes anything. Load does not
refuse a file it loaded before. **Save Session** still refuses what it refused
before.

Do not change a version number, resource hash, or runtime binding by hand;
changing metadata does not migrate its content.

## What a session preserves

A current session embeds its input resources, normalized settings, the last
committed render request and generated Result of each diagram mode, and
supported editor and comparison state. Replay does not depend on the original file path remaining valid.
Treat the session as sensitive when its embedded source data is sensitive.

Each diagram mode keeps its own feature label style, **Label Line** and
**Label Protection** (`--label_line`, `--label_protection`). A Session saved
before they existed loads with each mode's defaults: Linear Underline and Box,
Circular None and None.

An embedded resource may declare `checksum`, the SHA-256 digest of its bytes as
`sha256:<hex>` or as bare hex. gbdraw does not add one. The web app keeps it
when it saves the session again, and gbdraw rejects a resource whose bytes do
not match it.

Replay reads embedded BLAST outfmt 6/7 resources with the current [comparison
table rules](input-formats-and-tsv-schemas.md#comparison-and-numeric-tables).
A saved table with more than 12 columns replays from its first 12 columns. A
malformed saved table, or one whose IDs name the wrong Linear records, stops
replay with an error.

Both `.gbdraw-session.json` and lossless `.gbdraw-session.json.gz` are
accepted. The web app writes compressed sessions by default. The command line
writes `<output>.gbdraw-session.json` by default. `--session_output` selects an
explicit path and implies saving; a `.gz` suffix selects compression.

The web app keeps one Result for each diagram mode. **Save Session** writes
every Result and the editable state: the shown mode's Result at the top level
when it has one, and the other mode's in `otherModeResult`. **Load Session**
restores saved values instead of applying fresh browser defaults, and shows the
saved mode when that mode has a Result, otherwise the mode that has one. Generate before you save when the Result should reflect changed controls.
Saving before Generate keeps the newer draft alongside the earlier Result.
Loading and regenerating keep biological identities, labels, record
placement, comparison artifacts, and supported editor state. SVG bytes or text
metrics can still differ across gbdraw versions.

A Web Session also preserves supported manual legend, plot title, and Linear
scale positions in its saved Result. After loading, **Generate Diagram** carries
each matched item's offset into the new automatic layout for the same diagram;
it does not pin an absolute page position. If the source, selected region,
record set, or layout no longer matches a moved item, Generate keeps the saved
Result and asks for that item's position reset, **Reset Layout**, or matching
settings. Preview search placement, whether the search bar is hidden, and
Layout edit hints are not saved. Keep the Session when browser positioning must
be reproduced; a raw Python render request alone does not contain those manual
positions. Session version 40 files that unreleased builds of early August 2026
saved recorded moved items only in the drawing: Load shows them where they
were, and Generate lays them out automatically.

Fresh CLI sessions, and sessions that the Python API saves, omit `config`
because they have no independent Web draft. Web initializes their settings from
`renderRequest`, including the tables it records (`-t`, `-d`,
`--feature_visibility_table`, `--label_table`, `--label_whitelist`,
`--qualifier_priority`, and `--feature_override_table`), so the next
**Generate Diagram** keeps them. Web restores the original input files from
their bindings. A Session version 40 or 41 that the CLI wrote carries a
`config` derived from its options; Web and CLI replay read it as no draft, so
these Sessions load the same way. A Web draft in `config` must contain valid
`form` and `adv` objects; a partial draft is rejected. CLI replay preserves a
supplied Web draft.

Which records Web draws from a CLI session depends on the request:

- When the request draws every record of an input file, its record source names
  that file's resource, so feature-popup record rotation works right after
  Load.
- A Linear request, or a single Circular request, can draw only some records of
  a multi-record file (for example with `--record_id` or a records table). It
  also reads that file and selects each drawn record by record ID, or by `#n`
  when another record of the file has the same ID. The file stays one File that
  draws only those records.
- A Linear record, or the record of a single Circular request, can be drawn
  cropped or reverse-complemented. The session then stores its crop in source
  coordinates as `region` and its orientation as `region.reverseComplement` or
  `presentation.reverseComplement`, as Web Save does. Each Linear row takes its crop and orientation from the request, so
  **Generate Diagram** draws the record the CLI drew.
- A cropped or reversed record of a Circular batch or grid, the records of a
  file that a Circular batch or grid draws only partly, and every record of an
  older CLI session keep their own drawn copy. Rotating such a record needs
  **Generate Diagram** first.

A CLI binding uid such as `cli-seq-1` is only an initial value: each Linear
file takes the record identity of `renderRequest.records[].recordKey`. A
multi-record file, whose records are `record-1:1`, `record-1:2`, and so on,
becomes one Linear row per record with that recordKey and the request's record
selector (the `#n` of its recordKey when an older CLI session stored no
selector), so **Inherit saved comparison** finds every record it names.

A CLI comparison stays read-only in Web. **Inherit saved comparison** reuses it
after promoting an older committed request, such as a version 42 sidecar, to
the current schema. A Linear Session without a stored comparison plan loads
with **No comparison**. A CLI Session written with `-b` therefore offers
**Replace with current controls** only after a comparison is set up, and starts
no LOSAT run before that.
**Generate Diagram** before choosing **Inherit saved comparison**, **Replace
with current controls**, or **Clear comparison** stops with
`COMPARISON_CHOICE_REQUIRED` ("Choose how to handle the saved comparison in the
Comparison panel, then Generate again."). A saved comparison whose resource is
missing reports `COMPARISON_INPUT` instead.
After **Replace with current controls**, Generate checks the comparison controls
as it does for an editable comparison, so an incomplete pair reports its
`COMPARISON_INPUT` issue. Controls that set no comparison stop with
`COMPARISON_REPLACEMENT_EMPTY` ("Set up a comparison in the Comparison panel, or
choose Clear comparison, then Generate again.").

A CLI Linear Session written without `-b` or `--losat` commits only a disabled
protein pipeline (mode `none`, no pairs). It loads with **No comparison** and
without selecting LOSATP, so Generate draws no ribbons and starts no LOSAT run.
A CLI Session written with `--losat losatp` and an editable protein pipeline
loads with the adjacent LOSATP comparison it drew. It stores the protein
results the CLI computed: Similarity groups stores the comparison tables and
the `orthogroupResult`, and Collinear stores only the `collinearityResult` (no
tables). The first Generate after Load reuses them instead of running LOSAT
again. Pairwise stores nothing, so its first Generate still searches.

Without saved `ui.layoutPreferences`, the legend and plot-title positions come
from the committed `diagramOptions.output` for its mode and grouping, so a CLI
`--legend` survives loading and the first Generate.

The Web file inventory uses binding schema 2. An explicit composite `c_gb`
restores one editable GenBank File from ordered component resource bindings.
Each component keeps its own filename, MIME type, modification time and exact
bytes; duplicate payloads can share storage while every occurrence is kept.
Combining appends LF to each nonempty component that lacks a final LF.
Combining does not depend on the last committed `renderRequest`; replay still
uses that request. Python validates and preserves the draft binding without combining
its components for rendering.

Schema-1 ordinary bindings and File arrays remain supported. Existing sessions
without explicit bindings keep their request-derived source initialization;
original components cannot be recovered if their membership was never saved.
Schema 2 is accepted with sessions 41–42, 44, and 46–48. Python refuses unknown or malformed
bindings; the Web app drops each one and lists it in the **Load this file?** dialog
before the import replaces the current work. Older schema-1 readers reject new
schema-2 documents; changing the schema number does not convert them.

## Fonts

A Session and a render request store the font by name (the CSS font-family
list, `objects.text.font_family`), never the font file. The Session's preview
never embeds a font; only exported SVG does (see
[Fonts](output-formats-and-export.md#fonts)). For a custom font, the web app
notes this under **Font** and in the **Save Session** notice, and
`--save_session` logs it.

A Session that stores only the full resolved configuration, as the command
line writes it, stores its font list in `objects.text.font_family`. When that
is the packaged list, the web app loads it as the default **Font** (`null`),
so **Save Session** then writes no font override; any other list loads as the
chosen **Font**.

When the first family of a Session's font list is not available, the command
line and Python log a warning and draw the nearest font: the first available
family of the list, else the Liberation face of its generic kind (serif,
sans-serif, or monospace). A later family that is not available changes
nothing, because it is not drawn. Give
the font with `--font_file` to use it. The web app asks before it loads such a
Session, without starting the drawing Worker: **Load font file…**, **Use
installed font…** (Chrome and Edge), **Load with** the nearest font, or
**Cancel**. **Load with** keeps the requested name in **Font**, shown as not
available, draws with the nearest font, and saves the original name again.

## Linear file-level defaults

A Linear source card carries a default **Organism / strain** and **Subtitle /
title** that apply to every record read from that file. A record that leaves its
own field empty inherits the file default; a record that fills it overrides.

`renderRequest.records[].presentation` always carries the resolved text that is
drawn. `webFiles.linearRecordMetadata[]` records the inheritance separately, so
loading restores the same editable state:

| Field | Meaning |
|---|---|
| `fileDefinition`, `fileSubtitle` | the source card's file-level defaults |
| `recordDefinition`, `recordSubtitle` | the record's own value, empty when it inherits |

Both pairs are written only when the source has a file-level default. Sessions
written before these fields existed carry neither; those readers compare the
resolved text against the file default instead, which reads an override that
happens to repeat the default as inheritance. Request schema 8 adds no field for
these Web-only defaults,
because these fields describe Web editing state rather than the
render request: a reader that ignores them still renders identical output.

`webFiles.bindings.linearSeqs[].inferred_definition` stores the definition
inferred from the record that row selects. It is written when a GenBank File is
uploaded or a row selects another record, and a record's definition is its own
value, then the file default, then this value. Sessions saved before the field
existed, and Sessions written by the CLI, carry none and are not re-inferred when
loaded, so their records keep the definitions they were drawn with.

`modes.<mode>.config.recordsOff` lists the records the Web app does not draw:
a Linear card uid, or a Circular source selector such as `#3`. Omission or `[]`
means every record is drawn. The Web app writes it only when a record is off, and Session version
46 is unchanged. The settings of an off record stay in the Session. The render
request holds only the records that are drawn, so it has no new field. A CLI
`--session` replay draws only those records, and a CLI re-save keeps
`recordsOff` and the settings of the off records.

Session 44 replaces the two editable Linear visibility booleans with independent
selected modes: `linear_accession_visibility` and `linear_length_visibility`,
each set to `auto`, `show`, or `hide`. Auto resolves from the effective rendered
rows and is projected to the existing request-schema-7 booleans; the request and
Python render model did not gain a new field. A version-42 `true` becomes Show,
`false` becomes Hide, and a missing boolean becomes historical Show. If a
selected-mode field is already present, it takes precedence. Current writers do
not write the retired booleans.

In a Session 27–44, the plot-title text (`plot_title`), plot-title font size
(`plot_title_font_size`), the record-label default font size (`def_font_size`),
and the comparison thresholds of the mode that was not shown are in
`config.modeProfiles`; the shown mode's values are the flat `form` and `adv`
fields. Load gives each mode its own values. When a Session has no per-mode
entry, its saved flat value belongs to the shown mode only, and the other mode
starts with its defaults. A missing `config.linearRecordLayout` means **Arrange in
rows** is on. A missing `ui.linearTypographyLinked` means the scale and ruler
label font sizes are linked while they are equal. A Circular request supplies no
Linear display values, so a Session without a saved draft keeps Accession and
Length at Auto and Replicon off.

When several records share a Linear row, a label or subtitle that no record of
the row contradicts describes the whole row and is drawn once beside it. An
empty value does not contradict anything, so a records table may name a row on
its leading record alone. A value that differs from the one the leading record
carries is record-local, and every record of the row then draws that line above
itself, including the record that leads the row. A subtitle follows its label,
so it is never left beside the row on its own.

Older Linear sessions may contain automatically inferred replicon or organelle
names saved as subtitles. These values remain visible when the **Replicon** row's **Show** box is
off: loading does not guess their origin or remove matching text. Clear an
unwanted record subtitle to restore its file default; clear that default too if
no subtitle is wanted. New uploads leave automatic names to the **Replicon** **Show** box.
Loading preserves the saved preview. **Generate** applies the current settings
and definition alignment, so its placement can differ from an older preview.

## Saving settings before loading a source

**Save Session** also works before any biological source is loaded. It preserves
the full editable configuration, Circular and Linear mode profiles, supported
editor preferences, and auxiliary files such as colors, filters and qualifier
priorities. **Load Session** replaces the current Session with those settings,
clearing any previous sources and Result. A rejected Load restores the previous
work. Load a real source and Generate to apply the saved settings to a diagram.

This settings-only variant was introduced in session 42 with explicit
`renderRequest: null`,
empty `results`, and a null feature catalog. It has no committed render. Missing
requests, dangling resources, or biological inputs in either mode cannot select
this variant. Auxiliary files retain their ordinary resource bindings and bytes.

Python can load and materialize a settings-only Session. CLI replay,
`session_to_request()` and `render_session()` report that it has no biological
render request. Existing supported full Sessions remain readable. Current
settings-only writers emit session 48, while current readers also accept the
session-42, session-44, session-46, and session-47 forms. Readers whose maximum
version is 42, 44, 46, or 47 reject newly written session-48 files.

Older settings JSON without a `format` field (containing `form` or `adv`) still
uses the legacy configuration import. It does not need a render request. This
is distinct from a malformed `format: "gbdraw-session", version: 41` envelope.

## Replay boundaries

On the command line, `gbdraw render --session FILE` replays every drawing that
has a Result (the top-level set and `otherModeResult`), or the drawings named
with `--drawing`. The `circular` and `linear` subcommands replay their own
mode's drawing and fail when the Session has no Result of that mode. Output
prefix, format, session-output, and overwrite options may replace their saved
counterparts. Other diagram options are rejected because they would combine
persisted and new settings ambiguously. The [command-line
reference](command-line.md#sessions-and-output) gives the output names.

When you replay a Session 31 or later with `--save_session` or
`--session_output`, the CLI first brings the Session to the current version and
validates it, so a Session that cannot be written again renders nothing. It
then replaces each rendered drawing's Result set where it is, at the top level
or in `otherModeResult`, keeps the other mode's set and the saved `ui.mode`,
and preserves the editable draft's component bytes, order and File metadata.
The re-save also replaces what both sets share: `losatCache` holds the entries
the new render returns, `losatDerivedCache` keeps the Session's own entries
unchanged (a CLI render builds none, and a Session saved by a fresh CLI render
has an empty derived cache), and `proteinIdentityManifest` and `legacyArtifacts`
come from the new render (with both drawings, `legacyArtifacts` come from the
Linear one). A Circular render, one diagram or a batch, returns
the Session's own `proteinIdentityManifest` and `legacyArtifacts`, so the
protein entries that the re-save keeps still resolve. A Linear LOSATP render
tries only the pending legacy protein entries of the record pairs it
compares; an entry it cannot match stays pending, with a warning, so a later
Web Generate can still use it. Legacy derived evidence is dropped once the
Session has current derived entries, or when the render used the last pending
entry. If
the render migrates older protein IDs (only Sessions 44 and older have them),
gbdraw rewrites the protein references throughout the Session. A resource
whose bytes are unchanged
keeps its ID and file name; a resource whose bytes changed is stored under a new
ID, and resources that nothing names any more are dropped. Explicit Web
bindings, including null and empty lists, take precedence over historical
direct-source lists.

In Python, `render_session()` is the persisted-session entry point.
`load_session_document()` validates a document, and `materialize_session()`
exposes embedded resources only while its context is active.
`session_to_request()` converts a current materialized session to a typed
request. Rendering that request alone does not replay saved comparison
artifacts; use `render_session()` when those artifacts belong in the result.

`render_request()` accepts current typed requests, not historical session
envelopes. Public typed session conversion accepts full versions 31–33, 39–42, 44, and 46–48;
versions 27–30 have no `renderRequest`: the CLI replays them, and the Web app
loads them from their saved settings. A Linear replay that combines `--region`,
`--reverse_complement 1`, and `-b` or `--comparisons_table` now crops before it
reverse-complements, so it can draw another window than the one its stored
comparison was searched on; it prints a `WARNING:` line, and the comparison may
need to be run again. Request schema 9 keeps schema 7's
display values and schema 6's input cardinality, including selectorless `all`
inputs. Resolve a typed request
before encoding when it still contains deferred paths or collection-level transforms.

## Similarity alignment in requests and Sessions

For Linear requests, `renderRequest` schemas 8 and 9 store `recordTranslations`
and `similarityAlignment` inside `renderRequest.layout`. The nested alignment
plan is schema 2; this does not change the Session version or request schema.
Every translation has one stable `recordKey` and finite `x` and `y` values. An
active plan covers those displayed keys. It stores the exact reference and each
target's Select or Skip outcome and rationale. The record's presentation or
region setting decides the orientation used to project each anchor center.
Current readers reject partial, malformed, mismatched, or unsupported plans.
Schema 1 of the nested plan was never released and has no reader. There is no
Circular form or generic transform matrix. A Python
`SimilarityAlignmentReference` is not saved: writers store the plan it resolves
to.

Web **Align…** and **Review alignment options…** read the direct ortholog links
for the selected Similarity Group from the orthogroup result in the committed
`renderRequest`, not from the feature catalog. A link counts only when both
ends bind to current members of their own source records. Missing or unbound
links leave the usual Review rules in place instead of a guess. A Web LOSATP
Generate commits its typed orthogroup or Collinear result with the request, so
Align, record rotation, and Save and Load use the same links without running
LOSATP again.

Released request schemas 1, 2, 5, 6, and 7 remain readable. Their
`alignOrthogroupFeature` protein-setting string is read only when loading old
files. Current writers emit neither that field nor the old Session-only
`orthogroupState.selectedOrthogroupAlignmentFeature` copy. A released older
Session with saved feature-catalog and orthogroup identity metadata is
converted to the current schema-2 plan before saving. Malformed or unmappable
old values produce an actionable error without replacing the last successful
Result. Current requests reject group-only input. Loading a saved preview does
not start LOSATP; [Save and Load Sessions](web-app.md#save-and-load-sessions)
describes when Load starts the diagram engine.

A current Session round trip keeps exact feature and record identities,
Select/Skip rationale, record orientation, base translations, and the data
that **Reset alignment…** needs: the immediate pre-align positions and the
directions that Align changed. Ordinary Generate renders the saved plan without
guessing a new anchor. A stable reorder resolves by `recordKey` and biological
feature identity. Source replacement, crop, selector, and record drag clear the
plan with a visible reason. Manual Reverse keeps the plan and aligns the same
anchors on the next Generate.

**Reset alignment…** offers positions only (keeping current directions) or
positions plus the latest Align's actual direction changes. Combined Reset
restores the absolute pre-Align directions only for those records, including a
changed reference, and replaces later manual direction edits on those targets.
Either scope clears the plan and the saved Reset data; Undo restores them
before you can choose another scope. See
[Web Reset](web-app.md#reset-positions-or-directions).

A stale reference requires Reselect/Clear, and a stale target requires
Select/Skip. While a repair is pending or has failed, the last successful
Result stays. Undo/Redo restores the whole alignment. Preview-only guides,
candidate markers, and recommendation badges are never saved.

The Reset data is stored in `editorState.alignmentResetReceipt`, separately
from the direction-independent plan. It binds source bytes, record/selector/crop
identity and exact anchors. It stores only actual absolute direction deltas
plus any reference x correction, and never a review mode or Custom policy. Save
and fresh Load keep it. Style regeneration, stable reorder and ordinary Reverse
do not rewrite the pre-Align directions. Source/selector/crop invalidation and
record drag clear it with the plan.

An older Session without saved direction data keeps its position baseline,
and combined Reset is disabled. When the saved Reset data has an empty
direction list, the dialog instead says that the latest Align changed no
directions. Malformed or stale Reset data is rejected before it replaces the
current alignment; it is not silently dropped. The Reset data needs no format
change of its own: the plan stays schema 2.

In Web **Run Info**, **Source recipe** uses the original input filenames and
public CLI settings. Keep those original files and download any listed
generated helpers with **Download reproducibility files**. That bundle also
includes the session that **Exact replay** refers to, with its embedded
resources and saved analysis results. An unavailable Source recipe reports why
it cannot express the committed settings losslessly.

Both commands target the successful generation represented by Run Info, even
after Undo/Redo or edits to the current controls. Exact replay reconstructs that
committed generation; it does not apply later preview-only editor changes or an
ungenerated draft. **Save Session** is the route for preserving editable work
alongside the earlier Result. Exact replay does not promise byte-identical SVG
across gbdraw versions or font environments.

## Saved comparison results and cache reuse

Sessions can retain raw comparison rows, derived Similarity groups or
Collinear blocks, and the protein identity information needed to bind them to
features. Every derived runtime handle must resolve to the identity information
stored with that artifact. A missing or inconsistent identity fails instead of
falling back to a display label.

Comparison cache reuse requires the same sequence content, selected proteins,
record and feature bindings, query/subject direction, program, and meaningful
search arguments. Filenames and display labels do not define cache identity.
Only affected record pairs rerun when one valid cache key changes. A Session
stores each raw result once: rows that share a cache key, such as two ring
files with one sequence, share one entry named by the first row.

Pairwise hit limits, Similarity-group member limits, and Collinear block
settings are derived options. Changing one recomputes the affected derived
result while retaining eligible raw search rows. Current derived artifacts
record their upstream raw keys and requested settings. Collinear results
also record the effective `cds` or `locus` unit kinds produced by `auto` and
whether orthogroup inference was enabled. Sessions preserve each LOSATP mode's
hit-limit drafts. A saved Collinear pipeline without `collinearInferOrthogroups`
uses its historical inference behavior (ON); new Web state defaults to OFF.
In the Web, any derived result read from a Session loads unchanged and is
recomputed once from the raw search rows when the analysis runs again.
The Web, the CLI, and the Python API compute the cache key of a derived result
the same way, from every option that changes the result.

**Save Raw LOSAT TSV** resolves generated protein handles to stable readable
aliases. Uploaded comparison TSV is never rewritten. Export raw results to keep
a lasting record; the cache exists to avoid repeated work.

Release notes record when support changed;
this page documents current support.

## Retired inputs

New commands and Python calls reject the names and values below. A Session or
a canonical request schema 1–2 saved by an earlier release may still use them:
the CLI and the Python API rewrite them before they draw. The Web app reads
those a Session 31–33 request stores the same way, except a Circular slot
`spacing` (see below the table). A retired CLI flag exits with status 2 and
names its replacement. A retired Python field raises `TypeError`; there is
no alias. Field names stored inside a Session file do not change.

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

The config value `labels.stroke_width` (the leader line width) is removed:
leaders and underlines are 0.07 em of the label font size wide, at least 1 px.
A `config_overrides` entry or `--config_override` for it fails with the
unknown-path error. A Session or stored config that still holds it loads, and
its diagram draws the new width.

Long options of more than one word use underscores, with two aliases that still
work: `--annotation-table` for `--annotation_table`, and
`--gc_content_tick_interval` for `--gc_content_large_tick_interval`.

In the CLI and the Python API, a Circular slot whose spacing an earlier release
stored as a factor still draws. Set explicit pixel gaps (`inner_gap_px`,
`outer_gap_px`) before you save it in the current format; a factor cannot be
saved without loss. The Web app drops a factor spacing of a Session 31–33 at
Load and reads a spacing without a unit as pixels; it lists both in
**Load this file?**.

## Record rotation and feature placement

Session 41 and request schema 7 add requested record display and feature placement
intent together. Schema 6 keeps its original cardinality and row-inheritance
meaning; session 40 keeps its committed request and editable config as the sources of replay and draft.

Each schema-7 record has `display: {isCircular, startCoordinate}`. Both values
are nullable. A null start and an explicit source coordinate 1 remain distinct.
`diagramOptions.featurePlacements` contains sorted exact record/biological-feature
identities with a Main or mode-compatible lane-1 target. Auto removes an override.
The tolerance is `canvas.feature_overlap_tolerance_bp`, a non-negative integer
with default 0. Resolved lanes, pixel coordinates and display fragments are not
requested persistence fields.

Supported older requests have unset display, empty placements and tolerance 0.
Saving a schema-5/6 Web session promotes it without requiring Generate. Versions
27–30 still load without a `renderRequest` (CLI replay and Web Load), and
unknown or development-only versions remain rejected.

Editable Web rotation and placement drafts are saved in config, separately from
the last successful request and Result. An inactive rotation start stays in the
draft but is omitted from the effective request. Load restores both states;
Generate applies the draft. Source replacement invalidates the replaced source's
bindings, even when its filename or input-card UID is unchanged.
