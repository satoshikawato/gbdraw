[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | **FAQ** | [About](./ABOUT.md)

# Frequently asked questions

## Is there a web GUI? Do I need Streamlit?

Use [https://gbdraw.app/](https://gbdraw.app/) for the hosted app, or run `gbdraw gui` locally after installation. Streamlit is not required. Local GUI analysis runs on your machine; the interactive gallery examples are hosted at [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/).

## Why do my CLI and browser renders differ slightly?

Small differences in label placement and legend sizing are expected. The CLI uses kerning-aware font metrics, while the web UI uses browser text metrics.

## How do I hide the coordinate scale without hiding the genome axis?

Add `--hide_scale` in either Circular or Linear mode. In Linear mode, this
removes the bottom scale bar or ruler and any coordinate ticks and labels drawn
on record axes, while retaining each record's main axis line. In Circular
simple layouts, it removes the primary coordinate ticks and labels while
retaining the circular axis.

An explicit Circular track list is authoritative. If it contains an enabled
`ticks` slot, that slot remains visible even with `--hide_scale`; omit or
disable the slot instead. In the web app, clear **Show Coordinate Scale** for a
simple layout. While **Custom Track Slots** are active, use the **Ticks** slot
instead.

GC-content and Depth axes and ticks are separate controls. Linear definition
text such as **Length / Coord.** is also independent; use `--hide_length` to
remove that line.

## Can I use a GFF3 file by itself?

No. `gbdraw` requires both annotation and sequence data. When using GFF3 input, provide the matching FASTA file with `--fasta`.

```bash
gbdraw circular --gff my_genome.gff --fasta my_genome.fasta -o my_plot
```

## My labels overlap. What should I do?

Common fixes:

1. Reduce `--label_font_size`
2. Hide noisy labels with `--label_blacklist`
3. Keep only important labels with `--label_whitelist` regex patterns
4. Use the `--track_type middle` circular preset or reduce the number of displayed labels

See [Set feature colors and labels](./TUTORIALS/3_Advanced_Customization.md) for examples.

## How do I change the color of one specific gene?

Use a feature-specific color table with `-t`. This matches selected features by qualifier values and assigns a color and legend label.

See [Set feature colors and labels](./TUTORIALS/3_Advanced_Customization.md) and [Recipes](./RECIPES.md).

## Why does a web app color edit create a qualifier rule for some labels and `hash` rules for others?

When you choose **Apply to all label** or **Apply to all source label**, the web app first tries to represent the selected group as one exact feature-specific rule. For example, CDS features labeled `wsv360-like protein` can become a `product` rule with the pattern `^wsv360-like protein$`. The app does this only when every selected feature resolves to the same feature type, qualifier, and value; that rule matches exactly the selected group among the loaded features; and existing specific rules do not make precedence ambiguous. Regex metacharacters in the value are escaped before the pattern is anchored, so the exact-label selection does not accidentally become a broader regex.

If a single qualifier rule would be unsafe—for example, because the displayed label was manually edited, different features obtain that label from different qualifiers, or the rule would also match an unselected feature—the app keeps one `hash` rule per feature. `hash` is normally an exact feature selector rather than a GenBank qualifier. This fallback preserves the selected scope when the features have distinct generated identities. Both forms appear under **Specific Rules (-t)** and are sent to **Generate Diagram**, which keeps regenerated colors and legends consistent with the preview.

Exact duplicate records are a special case: identical feature type, record ID, coordinates, and strand produce the same generated `hash`. The editor can distinguish their rendered `_record_N` instances in the current preview, but **Generate Diagram** cannot reconstruct a one-instance-only color rule from identical biological identities. Use a distinguishing qualifier when available, or edit the duplicate instances together.

## How do I mark a coordinate range or a group of features?

Use `--annotation_table` and bind each `set_id` to an `annotations` custom track slot. Coordinate targets are 1-based and inclusive. Feature targets use the qualifiers already loaded from GenBank or GFF3, such as `locus_tag=ABC_001`.

The same table works in Circular and Linear mode. See [Region annotation tables](./TUTORIALS/5_Table_Driven_Inputs.md#7-region-annotation-tables).

## My comparative diagram has no ribbons. What is usually wrong?

The most common causes are:

1. The BLAST file is not in outfmt 6 or 7
2. The BLAST file order does not match the genome input order
3. Filtering thresholds such as `--evalue`, `--bitscore`, `--identity`, or `--alignment_length` are too strict

See [Draw genome comparison links from precomputed BLAST results](./TUTORIALS/2_Comparative_Genomics.md) for a working example.

## How do I draw several Linear records without comparing them in the web app?

Select **No comparison** above the Linear input rows, then click **Generate
Diagram**. The app skips LOSAT and does not use an uploaded BLAST TSV for that
render. Uploaded comparison files and custom raw-result filenames can remain in
the saved session as inactive drafts; use their separate reuse actions before
they can participate in a later comparison.

See the screenshot in [Web app comparison plans](./TUTORIALS/2_Comparative_Genomics.md#web-app-comparison-plans)
for the exact control to select.

Use **Run LOSAT** or **Upload BLAST TSV** for all positional adjacent gaps. In
upload mode, a gap without a file is deliberately skipped. Use **Selected
pairs** when one diagram needs a mixture of LOSAT edges, uploaded edges, and
omitted edges. Edit an entry under **Adjacent gaps**, or click **Add** under
**Selected pairs and retained drafts**, to create that plan. An included
uploaded edge must have an active file.

## Why did gbdraw rerun LOSATP after I loaded a session?

A protein-search cache hit requires the same amino-acid sequences, selected
proteins, record and feature bindings, direction, program, and meaningful search
arguments. A changed filename or display label does not invalidate the raw
search. Changed biological input or search settings do, and only the affected
record pair is rerun.

## Why does Save Raw LOSAT TSV not contain the internal `h_` IDs?

Those handles are session-internal references. For a generated protein result,
**Save Raw LOSAT TSV** resolves them to readable protein or feature aliases. If
any handle cannot be resolved safely, the download fails rather than exposing
an internal ID. Uploaded comparison TSV is never rewritten.

The full cache and export contract is documented in
[Session and request compatibility](./SESSION_COMPATIBILITY.md#saved-protein-comparison-results).

## Why is there less empty space between Linear comparison rows?

Automatic spacing now treats record bodies, comparison corridors, and definition text as separate X-aware constraints and uses the largest required clearance. A left-side definition block is not added to a plot-column comparison corridor when their horizontal ranges do not overlap, and `--comparison_height` is reserved only at boundaries crossed by a comparison. The value remains a minimum clear corridor; dense labels or tracks can still require more space.

## What if one record has no Depth TSV for a sample?

Keep the record's position in the repeated `--depth_track` group with a quoted empty argument:

```bash
gbdraw linear \
  --gbk record-a.gb record-b.gb \
  --depth_track record-a.depth.tsv '' \
  -o depth-partial \
  -f svg
```

Use `--depth_track '' record-b.depth.tsv` when only the second record has data. The empty argument means that the logical series is missing for that record. gbdraw does not substitute another file or draw zero coverage. In Linear mode the missing cell also reserves no vertical space, while the series identity remains stable for records that do contain data. Each group must contain at least one real file. See [Plot read depth and other numeric tracks](./TUTORIALS/6_Depth_Quantitative_Tracks.md#3-compare-depth-across-records) for a runnable example.

## Why is my circular BLAST similarity ring empty?

Check that the BLAST file is outfmt 6 or 7, the displayed circular record ID appears on the side selected by `--conservation_reference`, and the thresholds are not too strict. When BLAST was generated as `blastn -query comparison.fasta -subject reference.fasta`, use `--conservation_reference subject`.

These rings draw raw HSP spans; they do not infer evolutionary conservation. A BLAST row where the selected reference start is greater than the selected reference end is treated as reverse orientation, not as a hit crossing the circular origin. The current implementation does not infer binned or wraparound hits.

## Can pairwise comparison links be curved?

Yes. In linear mode, `--pairwise_match_style ribbon` draws straight filled ribbons by default. Use `--pairwise_match_style curve` to bend the same match spans; curved links can be easier to distinguish in a dense comparison diagram.

## Can I use gene names instead of product descriptions for labels?

Yes. Use `--qualifier_priority` to prefer `gene`, `locus_tag`, or other qualifiers.

```tsv
CDS	gene
```

```bash
gbdraw circular --gbk genome.gb --labels --qualifier_priority priority.tsv -o output -f svg
```

## How do I make the GC content track smoother?

Increase the window and step sizes:

```bash
gbdraw circular --gbk genome.gb --window 10000 --step 1000 -o output -f svg
```

The montage uses one *E. coli* record and changes only `--window`/`--step`: 100000/10000, 10000/1000, and 1000/100 from left to right. This wider range makes the smoothing tradeoff visible without removing the feature, metadata, legend, GC content, or GC skew context.

![Three E. coli diagrams showing coarse, medium, and fine GC-content and GC-skew windows](../examples/window_step_comparison.png)

## Can I plot AT instead of GC?

Yes. Use `--nt AT`.

```bash
gbdraw circular --gbk genome.gb --nt AT -o output -f svg
```

The 12-panel comparison keeps the *E. coli* record, feature tracks, metadata, and legend fixed. In reading order it shows GC, CG, AG, GA, CT, TC, TG, GT, CA, AC, AT, and TA; only `--nt` changes between panels. Reversing a pair keeps its content track unchanged while reversing the sign of its skew.

![Twelve E. coli diagrams comparing dinucleotide content and skew for GC, CG, AG, GA, CT, TC, TG, GT, CA, AC, AT, and TA](../examples/skew_comparison.png)

## Why does SVG export work but PNG/PDF/EPS/PS export fail?

Non-SVG export requires CairoSVG. Install the optional export dependency and, if needed on your platform, the system Cairo/Pango libraries.

## Are there known visualization limitations?

- Trans-introns are not currently visualized.
- Mixed-format text such as `<i>Ca.</i> Tyloplasma litorale` does not reliably survive SVG-to-PNG/PDF/EPS/PS conversion. Use SVG if you need exact mixed formatting.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | [Gallery](./GALLERY.md) | **FAQ** | [About](./ABOUT.md)
