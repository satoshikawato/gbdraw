[Documentation home](../DOCS.md) | [Web app Tutorial](../TUTORIALS/GUI/README.md) | [Web app how-to guides](../HOW_TO/GUI/README.md)

# Web app reference

The hosted application at [gbdraw.app](https://gbdraw.app/) and local `gbdraw gui` expose the same single-page interface. The hosted page loads static assets from the site; parsing, comparison search, rendering, session handling, and export run in the browser.

## Modes and inputs

| Mode | Primary input | Result ownership |
|---|---|---|
| Circular | One GenBank/GBFF file or one matched GFF3 + FASTA pair | One result, a multi-record canvas, or separate results according to grouping |
| Linear | GenBank rows or matched GFF3 + FASTA rows | One ordered Linear result with an independent comparison plan |

A GenBank file may contain several records. Circular input from several separate source files belongs in the CLI or Python API. Linear input cards own their selector, region, reverse-complement flag, record label, subtitle, and grid placement.

## Main workflow

1. Select **Circular** or **Linear**.
2. Add biological input files and resolve any validation message.
3. Choose layout, track, comparison, label, and output settings.
4. Select **Generate Diagram**.
5. Inspect **Result Preview**, then use a format button or save a session.

Settings are not committed to a result until generation succeeds. Preview edits and saved draft controls may exist beyond the last committed render, so a session stores both the committed request and supported editable state.

## Browser defaults

| Setting | Circular | Linear |
|---|---|---|
| Separate strands | On | On |
| Legend | Left | Bottom |
| Feature placement | Tuckin preset | Features on axis |
| Comparison | Not applicable until similarity rings are configured | No comparison until a plan is selected |

A loaded session restores its saved values instead of applying these defaults. Custom track stacks keep their draft slots while inactive; **Reset** is the action that rebuilds them from simple controls.

## Comparison surfaces

The browser can run LOSATN, TLOSATX, and LOSATP. LOSATP supports Pairwise protein links, Similarity groups, and Collinear blocks. Uploaded BLAST-compatible tables can provide prepared Linear links or Circular similarity-ring evidence. The Linear comparison plan is independent of row placement and supports no comparison, adjacent comparisons, or explicit selected edges.

See [Comparison programs, thresholds, and result semantics](comparison-programs-thresholds-and-results.md) for method boundaries and retained-result meanings.

## Preview and export

The preview exposes semantic records, tracks, features, comparisons, definitions, ticks, legends, and annotations. Interactive output may add search, feature and match popups, group inspection, sequence downloads, and embedded metadata. Available export actions are described in [Output formats and export](output-formats-and-export.md).

## Accessibility and automation

Primary Tutorial controls have stable accessible names: **Circular**, **Linear**, **GenBank/DDBJ File**, **Output Prefix**, **Species**, **Track Preset**, **Separate Strands**, **Hide GC Content**, **Hide GC Skew**, **Label Mode**, **Legend Position**, **Generate Diagram**, **Result Preview**, and **SVG**. Mode buttons expose pressed state, file controls are labelled, and preview output is a named region.

File-size and runtime limits depend on browser memory, feature density, comparison topology, and search settings. The application validates known structural limits but does not promise one universal maximum genome size. See [Browser privacy, offline execution, and performance](../EXPLANATION/browser-privacy-offline-execution-and-performance.md).
