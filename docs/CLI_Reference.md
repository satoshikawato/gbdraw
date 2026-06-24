[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | **CLI Reference** | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)

# Command-Line Reference

This reference mirrors the current command help from `python -m gbdraw.cli`. Use it as the source of truth for available options and defaults.

## Main Command

```text
gbdraw v. 0.9.2: A diagram generator for small genomes

Usage:
  gbdraw <subcommand> [options]

Subcommands:
  circular  Generate a circular genome diagram
  linear    Generate a linear genome diagram
  gui       Launch the local web UI in your browser

For each subcommand, you can get additional help by running:
  gbdraw <subcommand> --help

Examples:
  gbdraw circular --gbk input.gb
  gbdraw circular --gff input.gff --fasta input.fna
  gbdraw linear --gbk input.gb
  gbdraw linear --gff input.gff --fasta input.fna
  gbdraw linear --gbk input1.gb input2.gb input3.gb -b input1_input2.blast.outfmt7.txt input2_input3.blast.outfmt7.txt
  gbdraw linear --gff input1.gff input2.gff input3.gff --fasta input1.fna input2.fna input3.fna -b input1_input2.blast.outfmt7.txt input2_input3.blast.outfmt7.txt
  gbdraw gui

Options (examples):
  --gbk                Input GenBank file(s)
  --gff                Input GFF# file(s) (rquires --fasta; mutually exclusive with --gbk)
  --fasta              Input FASTA file(s) (required with --gff; mutually exclusive with --gbk)
  -o, --output         Output file prefix (optional)
  -b, --blast          BLAST result file in tab-separated format (-outfmt 6 or 7) (optional; implemented for linear mode only)
  --conservation_blast BLAST result file(s) for circular conservation rings (-outfmt 6 or 7)

Additional Information:
  - For full documentation, visit: https://github.com/satoshikawato/gbdraw/
  - For issues and source code, visit the GitHub repository: https://github.com/satoshikawato/gbdraw/
  - For support, contact: kawato[at]kaiyodai.ac.jp
```

`gbdraw gui` starts a local HTTP server on a free port and opens the web UI in your browser. Streamlit is not required.

## Circular Mode

```text
usage: cli.py [-h] [--gbk [GBK_FILE ...]] [--gff [GFF3_FILE ...]]
              [--fasta [FASTA_FILE ...]] [-o OUTPUT] [-p PALETTE] [-t TABLE]
              [-d DEFAULT_COLORS] [-n NT] [-w WINDOW] [-s STEP]
              [--species SPECIES] [--strain STRAIN] [-k FEATURES]
              [--feature_shape TYPE=SHAPE]
              [--block_stroke_color BLOCK_STROKE_COLOR]
              [--block_stroke_width BLOCK_STROKE_WIDTH]
              [--axis_stroke_color AXIS_STROKE_COLOR]
              [--axis_stroke_width AXIS_STROKE_WIDTH]
              [--line_stroke_color LINE_STROKE_COLOR]
              [--line_stroke_width LINE_STROKE_WIDTH]
              [--definition_font_size DEFINITION_FONT_SIZE]
              [--plot_title PLOT_TITLE]
              [--plot_title_font_size PLOT_TITLE_FONT_SIZE]
              [--keep_full_definition_with_plot_title]
              [--center_reserved_radius CENTER_RESERVED_RADIUS]
              [--label_font_size LABEL_FONT_SIZE] [-f FORMAT] [--suppress_gc]
              [--suppress_skew]
              [--conservation_blast BLAST [BLAST ...]]
              [--conservation_reference {query,subject,auto}]
              [--conservation_labels LABEL [LABEL ...]]
              [--conservation_ring_width CONSERVATION_RING_WIDTH]
              [--conservation_ring_gap CONSERVATION_RING_GAP]
              [--evalue EVALUE] [--bitscore BITSCORE] [--identity IDENTITY]
              [--alignment_length ALIGNMENT_LENGTH]
              [-l LEGEND] [--multi_record_canvas]
              [--multi_record_size_mode {auto,linear,equal,sqrt}]
              [--multi_record_min_radius_ratio MULTI_RECORD_MIN_RADIUS_RATIO]
              [--multi_record_column_gap_ratio MULTI_RECORD_COLUMN_GAP_RATIO]
              [--multi_record_row_gap_ratio MULTI_RECORD_ROW_GAP_RATIO]
              [--multi_record_position MULTI_RECORD_POSITION]
              [--plot_title_position {none,top,bottom}] [--separate_strands]
              [--track_type TRACK_TYPE] [--resolve_overlaps]
              [--labels [{none,out,both}]]
              [--label_rendering {auto,embedded_only,external_only}]
              [--label_whitelist LABEL_WHITELIST |
              --label_blacklist LABEL_BLACKLIST]
              [--qualifier_priority QUALIFIER_PRIORITY]
              [--label_table LABEL_TABLE] [--feature_table FEATURE_TABLE]
              [--outer_label_x_radius_offset OUTER_LABEL_X_RADIUS_OFFSET]
              [--outer_label_y_radius_offset OUTER_LABEL_Y_RADIUS_OFFSET]
              [--inner_label_x_radius_offset INNER_LABEL_X_RADIUS_OFFSET]
              [--inner_label_y_radius_offset INNER_LABEL_Y_RADIUS_OFFSET]
              [--scale_interval SCALE_INTERVAL]
              [--feature_width FEATURE_WIDTH]
              [--circular_track_order CIRCULAR_TRACK_ORDER]
              [--circular_track_slot CIRCULAR_TRACK_SLOT]
              [--gc_content_width GC_CONTENT_WIDTH]
              [--gc_content_radius GC_CONTENT_RADIUS]
              [--gc_content_mode {deviation,percent}]
              [--gc_content_min_percent GC_CONTENT_MIN_PERCENT]
              [--gc_content_max_percent GC_CONTENT_MAX_PERCENT]
              [--gc_content_tick_interval GC_CONTENT_TICK_INTERVAL]
              [--gc_content_large_tick_interval GC_CONTENT_LARGE_TICK_INTERVAL]
              [--gc_content_small_tick_interval GC_CONTENT_SMALL_TICK_INTERVAL]
              [--gc_content_tick_font_size GC_CONTENT_TICK_FONT_SIZE]
              [--show_gc_content_axis] [--hide_gc_content_axis]
              [--show_gc_content_ticks] [--hide_gc_content_ticks]
              [--gc_skew_width GC_SKEW_WIDTH]
              [--gc_skew_radius GC_SKEW_RADIUS]
              [--legend_box_size LEGEND_BOX_SIZE]
              [--legend_font_size LEGEND_FONT_SIZE]

Generate genome diagrams in PNG/PDF/SVG/PS/EPS. By default, diagrams for
multiple entries are saved separately. Use --multi_record_canvas to place
multiple records on one grid canvas.

options:
  -h, --help            show this help message and exit
  --gbk [GBK_FILE ...]  Genbank/DDBJ flatfile
  --gff [GFF3_FILE ...]
                        GFF3 file (instead of --gbk; --fasta is required)
  --fasta [FASTA_FILE ...]
                        FASTA file (required with --gff)
  -o, --output OUTPUT   output file prefix (default: accession number of the
                        sequence)
  -p, --palette PALETTE
                        Palette name (default: default)
  -t, --table TABLE     color table (optional)
  -d, --default_colors DEFAULT_COLORS
                        TSV file that overrides the color palette (optional)
  -n, --nt NT           dinucleotide (default: GC).
  -w, --window WINDOW   window size (optional; default: 1kb for genomes < 1Mb,
                        10kb for genomes <10Mb, 100kb for genomes >=10Mb)
  -s, --step STEP       step size (optional; default: 100 bp for genomes <
                        1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb)
  --species SPECIES     Species name (optional; e.g. "<i>Escherichia
                        coli</i>", "<i>Ca.</i> Hepatoplasma crinochetorum")
  --strain STRAIN       Strain/isolate name (optional; e.g. "K-12", "Av")
  -k, --features FEATURES
                        Comma-separated list of feature keys to draw (default:
                        CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)
  --feature_shape TYPE=SHAPE
                        Feature shape override (repeatable): TYPE=SHAPE where
                        SHAPE is arrow or rectangle.
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "gray")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (optional; float; default: 2 pt for
                        genomes <= 50 kb, 0 pt for genomes >= 50 kb)
  --axis_stroke_color AXIS_STROKE_COLOR
                        Axis stroke color (str; default: "gray")
  --axis_stroke_width AXIS_STROKE_WIDTH
                        Axis stroke width (optional; float; default: 3 pt for
                        genomes <= 50 kb, 1 pt for genomes >= 50 kb)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (str; default: "gray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (optional; float; default: 5 pt for
                        genomes <= 50 kb, 1 pt for genomes >= 50 kb)
  --definition_font_size DEFINITION_FONT_SIZE
                        Definition font size (optional; default: 18)
  --plot_title PLOT_TITLE
                        Circular plot title shown when plot_title_position is
                        top or bottom (optional; defaults to species +
                        strain).
  --plot_title_font_size PLOT_TITLE_FONT_SIZE
                        Plot title font size for circular top/bottom title
                        layout (optional; default: 32).
  --keep_full_definition_with_plot_title
                        Keep the full centered record definition when a
                        circular plot title is shown (default: False).
  --center_reserved_radius CENTER_RESERVED_RADIUS
                        Override the centered definition reserved radius for
                        circular track packing (in px; must be >= 0).
  --label_font_size LABEL_FONT_SIZE
                        Label font size (optional; default: 14 (pt) for
                        genomes <= 50 kb, 8 for genomes >= 50 kb)
  -f, --format FORMAT   Comma-separated list of output file formats (svg,
                        interactive-svg, png, pdf, eps, ps; default: svg).
                        PNG/PDF/EPS/PS require CairoSVG.
  --suppress_gc         Suppress GC content track (default: False).
  --suppress_skew       Suppress GC skew track (default: False).
  --conservation_blast BLAST [BLAST ...]
                        Precomputed BLAST outfmt 6/7 file(s) for circular
                        conservation rings.
  --conservation_reference {query,subject,auto}
                        BLAST side containing displayed circular reference
                        coordinates. Use "subject" for BLAST generated as
                        comparison FASTA query against displayed genome
                        subject. Default: "auto".
  --conservation_labels LABEL [LABEL ...]
                        Labels for conservation rings, aligned by logical
                        source index.
  --conservation_ring_width CONSERVATION_RING_WIDTH
                        Conservation ring width for circular mode (in px; must
                        be > 0).
  --conservation_ring_gap CONSERVATION_RING_GAP
                        Conservation ring gap for circular mode (in px; must
                        be > 0).
  --evalue EVALUE       Maximum BLAST e-value retained for conservation rings
                        (default: 1e-5).
  --bitscore BITSCORE   Minimum BLAST bitscore retained for conservation rings
                        (default: 50).
  --identity IDENTITY   Minimum BLAST identity percentage retained for
                        conservation rings (default: 70).
  --alignment_length ALIGNMENT_LENGTH
                        Minimum BLAST alignment length retained for
                        conservation rings (default: 0).
  -l, --legend LEGEND   Legend position (default: "right"; "left", "right",
                        "top", "bottom", "upper_left", "upper_right",
                        "lower_left", "lower_right", "none")
  --multi_record_canvas
                        Place multiple records on one shared canvas using
                        automatic grid layout (default: False).
  --multi_record_size_mode {auto,linear,equal,sqrt}
                        Size mode for multi-record circular canvas ("auto",
                        "linear", "equal"; "sqrt" is accepted as an alias of
                        "auto"; default: "auto").
  --multi_record_min_radius_ratio MULTI_RECORD_MIN_RADIUS_RATIO
                        Minimum radius ratio for multi-record scaling (0 <
                        ratio <= 1; default: 0.55).
  --multi_record_column_gap_ratio MULTI_RECORD_COLUMN_GAP_RATIO
                        Additional horizontal gap ratio between visible
                        content bounds in each multi-record row (>= 0;
                        default: 0.10).
  --multi_record_row_gap_ratio MULTI_RECORD_ROW_GAP_RATIO
                        Additional gap ratio between multi-record row content
                        bounds (>= 0; default: 0.05).
  --multi_record_position MULTI_RECORD_POSITION
                        Record placement for multi-record canvas (repeatable):
                        <selector>@<row> where selector is #index or record_id
                        and row starts at 1.
  --plot_title_position {none,top,bottom}
                        Plot title position in circular mode ("none", "top",
                        "bottom"; default: "none").
  --separate_strands    Separate strands (default: False).
  --track_type TRACK_TYPE
                        Circular preset for legacy/simple layout. Choices:
                        "tuckin", "middle", "spreadout". Ignored when
                        explicit --circular_track_slot layouts are supplied.
  --resolve_overlaps    Resolve overlapping features by placing them on
                        separate tracks (default: False). Useful for plasmid
                        visualization.
  --labels [{none,out,both}]
                        Label placement mode: no argument or "out" (outside),
                        "both" (outside+inside), or "none" (hidden). Default:
                        "none".
  --label_rendering {auto,embedded_only,external_only}
                        Label rendering policy. "auto" embeds labels that fit
                        and routes the rest externally; "embedded_only" drops
                        external labels; "external_only" forces labels outside
                        feature bodies. Default: "auto".
  --label_whitelist LABEL_WHITELIST
                        Path to a TSV file for label whitelisting by regex
                        pattern (optional); mutually exclusive with
                        --label_blacklist
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label
                        blacklisting (optional); mutually exclusive with
                        --label_whitelist
  --qualifier_priority QUALIFIER_PRIORITY
                        Path to a TSV file defining qualifier priority for
                        labels (optional)
  --label_table LABEL_TABLE
                        Path to a TSV file defining post-filter label text
                        overrides (optional)
  --feature_table FEATURE_TABLE
                        Path to a TSV file defining per-feature visibility
                        overrides (optional)
  --outer_label_x_radius_offset OUTER_LABEL_X_RADIUS_OFFSET
                        Outer label x-radius offset factor (float; default
                        from config)
  --outer_label_y_radius_offset OUTER_LABEL_Y_RADIUS_OFFSET
                        Outer label y-radius offset factor (float; default
                        from config)
  --inner_label_x_radius_offset INNER_LABEL_X_RADIUS_OFFSET
                        Inner label x-radius offset factor (float; default
                        from config)
  --inner_label_y_radius_offset INNER_LABEL_Y_RADIUS_OFFSET
                        Inner label y-radius offset factor (float; default
                        from config)
  --scale_interval SCALE_INTERVAL
                        Manual scale interval for circular mode (in bp).
                        Overrides automatic calculation.
  --feature_width FEATURE_WIDTH
                        Feature track width for circular mode (in px; must be
                        > 0).
  --circular_track_order CIRCULAR_TRACK_ORDER
                        Comma-separated circular slot order. Omitted slot
                        geometry inherits the selected --track_type preset for
                        each record; explicit slot fields override it.
  --circular_track_slot CIRCULAR_TRACK_SLOT
                        Circular track slot spec:
                        <slot_id>:<renderer>@key=value,key=value. Can be
                        repeated. Use r=<radius>, w=<width>,
                        inner_gap_px=<px>, and outer_gap_px=<px>. The legacy
                        spacing=<scalar> field is still accepted as a
                        compatibility alias for both gaps. ri/ro/gap are
                        obsolete. side and z are slot fields. If r, w, gap,
                        side, or standard renderer params are
                        omitted for built-in slots, they inherit the active
                        --track_type preset at render time. Inside numeric/depth
                        slots with no explicit r or w auto-compress when needed
                        and never move outside automatically. z only controls
                        SVG layering.
  --gc_content_width GC_CONTENT_WIDTH
                        GC content track width for circular mode (in px; must
                        be > 0).
  --gc_content_radius GC_CONTENT_RADIUS
                        GC content track center radius for circular mode (as a
                        ratio of base radius; must be > 0).
  --gc_content_mode {deviation,percent}
                        GC content display mode. deviation keeps the existing
                        mean-centered track. percent draws absolute GC
                        percentage as a baseline area track.
  --gc_content_min_percent GC_CONTENT_MIN_PERCENT
                        Minimum GC percent for percent-mode clipping/axis.
  --gc_content_max_percent GC_CONTENT_MAX_PERCENT
                        Maximum GC percent for percent-mode clipping/axis.
  --gc_content_tick_interval GC_CONTENT_TICK_INTERVAL
                        GC content percent-mode large tick interval; alias for
                        --gc_content_large_tick_interval.
  --gc_content_large_tick_interval GC_CONTENT_LARGE_TICK_INTERVAL
                        GC content percent-mode large tick interval.
  --gc_content_small_tick_interval GC_CONTENT_SMALL_TICK_INTERVAL
                        GC content percent-mode small tick interval.
  --gc_content_tick_font_size GC_CONTENT_TICK_FONT_SIZE
                        GC content percent-mode tick label font size.
  --show_gc_content_axis, --hide_gc_content_axis
                        Show or hide the GC content percent-mode axis.
  --show_gc_content_ticks, --hide_gc_content_ticks
                        Show or hide GC content percent-mode ticks and labels.
  --gc_skew_width GC_SKEW_WIDTH
                        GC skew track width for circular mode (in px; must be
                        > 0).
  --gc_skew_radius GC_SKEW_RADIUS
                        GC skew track center radius for circular mode (as a
                        ratio of base radius; must be > 0).
  --depth DEPTH         Depth TSV file in samtools depth format. Implies
                        --show_depth.
  --depth_track DEPTH [DEPTH ...]
                        Repeatable logical depth track. In circular mode,
                        provide one file for a single record, or one file per
                        record when using --multi_record_canvas.
  --depth_track_label LABEL [LABEL ...]
                        Depth track label(s). Provide one label or one per
                        --depth_track.
  --depth_track_color COLOR [COLOR ...]
                        Depth track fill color(s). Provide one color or one
                        per --depth_track.
  --show_depth          Show depth coverage track. Required only when no depth
                        file option is supplied.
  --depth_width DEPTH_WIDTH
                        Depth track width for circular mode (in px; must be
                        > 0).
  --legend_box_size LEGEND_BOX_SIZE
                        Legend box size (optional; float; default: 24 (pixels,
                        96 dpi) for genomes <= 50 kb, 20 for genomes >= 50
                        kb).
  --legend_font_size LEGEND_FONT_SIZE
                        Legend font size (optional; float; default: 20 (pt)
                        for genomes <= 50 kb, 16 for genomes >= 50 kb).
```

Circular conservation rings use one ring per `--conservation_blast` source and a shared identity gradient legend. BLAST tables must be outfmt 6 or 7. Coordinates on the selected reference side are normalized from BLAST 1-based inclusive coordinates to drawing spans; `start > end` marks reverse orientation and is not interpreted as a circular-origin-spanning hit. The CLI does not run LOSAT for conservation rings, so provide precomputed BLAST output.

Depth tracks can be supplied with the legacy `--depth` option or the repeatable
`--depth_track` option. `--depth` keeps the single-track SVG IDs `depth` and
`depth_axis`. Multiple `--depth_track` groups render as `depth_1`,
`depth_2`, and so on. Each `--depth_track` group is one logical track; provide
one file to reuse it for every record, or one file per record.

## Linear Mode

```text
usage: cli.py [-h] [--gbk [GBK_FILE ...]] [--gff [GFF3_FILE ...]]
              [--fasta [FASTA_FILE ...]] [-b [BLAST ...]] [-t TABLE]
              [--losatp_bin LOSATP_BIN]
              [--ncbi_blastp_bin NCBI_BLASTP_BIN]
              [--losatp_threads LOSATP_THREADS]
              [--protein_blastp_mode {none,pairwise,orthogroup,collinear}]
              [--collinear_min_anchors COLLINEAR_MIN_ANCHORS]
              [--collinear_max_unit_gap COLLINEAR_MAX_UNIT_GAP]
              [--collinear_color_mode {average_identity,orientation,orientation_identity}]
              [-p PALETTE] [-d DEFAULT_COLORS] [-o OUTPUT] [-n NT] [-w WINDOW]
              [-s STEP] [--separate_strands] [--show_gc]
              [--gc_content_mode {deviation,percent}]
              [--gc_content_min_percent GC_CONTENT_MIN_PERCENT]
              [--gc_content_max_percent GC_CONTENT_MAX_PERCENT]
              [--gc_content_tick_interval GC_CONTENT_TICK_INTERVAL]
              [--gc_content_large_tick_interval GC_CONTENT_LARGE_TICK_INTERVAL]
              [--gc_content_small_tick_interval GC_CONTENT_SMALL_TICK_INTERVAL]
              [--gc_content_tick_font_size GC_CONTENT_TICK_FONT_SIZE]
              [--show_gc_content_axis] [--hide_gc_content_axis]
              [--show_gc_content_ticks] [--hide_gc_content_ticks]
              [--show_skew]
              [--align_center] [--evalue EVALUE] [--bitscore BITSCORE]
              [--identity IDENTITY] [--alignment_length ALIGNMENT_LENGTH]
              [--pairwise_match_style {ribbon,curve}]
              [-k FEATURES] [--feature_shape TYPE=SHAPE]
              [--block_stroke_color BLOCK_STROKE_COLOR]
              [--block_stroke_width BLOCK_STROKE_WIDTH]
              [--axis_stroke_color AXIS_STROKE_COLOR]
              [--axis_stroke_width AXIS_STROKE_WIDTH]
              [--line_stroke_color LINE_STROKE_COLOR]
              [--line_stroke_width LINE_STROKE_WIDTH]
              [--definition_font_size DEFINITION_FONT_SIZE]
              [--plot_title PLOT_TITLE]
              [--plot_title_position {center,top,bottom}]
              [--plot_title_font_size PLOT_TITLE_FONT_SIZE]
              [--record_label RECORD_LABEL]
              [--label_font_size LABEL_FONT_SIZE]
              [--label_placement {auto,above_feature}]
              [--label_rendering {auto,embedded_only,external_only}]
              [--label_rotation LABEL_ROTATION]
              [--track_layout {above,middle,below}] [--track_axis_gap AUTO|PX]
              [--linear_track_order LINEAR_TRACK_ORDER]
              [--linear_track_slot SLOT]
              [--linear_track_axis_index LINEAR_TRACK_AXIS_INDEX]
              [--ruler_on_axis] [-f FORMAT] [-l LEGEND]
              [--show_labels [{all,first,none}]] [--resolve_overlaps]
              [--label_whitelist LABEL_WHITELIST |
              --label_blacklist LABEL_BLACKLIST]
              [--qualifier_priority QUALIFIER_PRIORITY]
              [--label_table LABEL_TABLE] [--feature_table FEATURE_TABLE]
              [--feature_height FEATURE_HEIGHT] [--gc_height GC_HEIGHT]
              [--comparison_height COMPARISON_HEIGHT]
              [--scale_style {bar,ruler}]
              [--scale_stroke_color SCALE_STROKE_COLOR]
              [--scale_stroke_width SCALE_STROKE_WIDTH]
              [--scale_font_size SCALE_FONT_SIZE]
              [--ruler_label_font_size RULER_LABEL_FONT_SIZE]
              [--ruler_label_color RULER_LABEL_COLOR]
              [--scale_interval SCALE_INTERVAL]
              [--legend_box_size LEGEND_BOX_SIZE]
              [--legend_font_size LEGEND_FONT_SIZE] [--normalize_length]
              [--region REGION] [--record_id RECORD_ID]
              [--reverse_complement REVERSE_COMPLEMENT]

Generate plot in PNG/PDF/SVG/PS/EPS.

options:
  -h, --help            show this help message and exit
  --gbk [GBK_FILE ...]  Genbank/DDBJ flatfile
  --gff [GFF3_FILE ...]
                        GFF3 file (instead of --gbk; --fasta is required)
  --fasta [FASTA_FILE ...]
                        FASTA file (required with --gff)
  -b, --blast [BLAST ...]
                        input BLAST result file in tab-separated format
                        (-outfmt 6 or 7) (optional)
  --losatp_bin LOSATP_BIN
                        Native LOSAT executable for --protein_blastp_mode
                        pairwise/orthogroup/collinear (default: losat).
  --ncbi_blastp_bin NCBI_BLASTP_BIN
                        NCBI BLAST+ blastp executable for
                        --protein_blastp_mode pairwise/orthogroup/collinear
                        (default: use automatic runtime resolution).
  --losatp_threads LOSATP_THREADS
                        Threads passed to the selected protein blastp runtime
                        for
                        --protein_blastp_mode pairwise/orthogroup/collinear
                        (default: runtime default).
  --protein_blastp_mode {none,pairwise,orthogroup,collinear}
                        Protein blastp comparison mode: none, pairwise
                        adjacent ribbons, all-record Orthogroups, or
                        Collinear blocks (default: none).
  --collinear_min_anchors COLLINEAR_MIN_ANCHORS
                        Minimum anchors/genes required for a rendered
                        Collinear block. The default 1 includes singleton
                        links.
  --collinear_max_unit_gap, --collinear_max_gene_gap COLLINEAR_MAX_UNIT_GAP
                        Maximum unit gap between neighboring collinear anchors
                        (default: 0).
  --collinear_color_mode {average_identity,orientation,orientation_identity}
                        Collinear ribbon color mode: average_identity,
                        orientation, or orientation_identity (default:
                        orientation).
  -t, --table TABLE     color table (optional)
  -p, --palette PALETTE
                        Palette name (default: default)
  -d, --default_colors DEFAULT_COLORS
                        TSV file that overrides the color palette (optional)
  -o, --output OUTPUT   output file prefix (default: out)
  -n, --nt NT           dinucleotide skew (default: GC).
  -w, --window WINDOW   window size (optional; default: 1kb for genomes < 1Mb,
                        10kb for genomes <10Mb, 100kb for genomes >=10Mb)
  -s, --step STEP       step size (optional; default: 100 bp for genomes <
                        1Mb, 1kb for genomes <10Mb, 10kb for genomes >=10Mb)
  --separate_strands    separate forward and reverse strands (default: False).
                        Features of undefined strands are shown on the forward
                        strand.
  --show_gc             plot GC content below genome (default: False).
  --gc_content_mode {deviation,percent}
                        GC content display mode. deviation keeps the existing
                        mean-centered track. percent draws absolute GC
                        percentage as a baseline area track.
  --gc_content_min_percent GC_CONTENT_MIN_PERCENT
                        Minimum GC percent for percent-mode clipping/axis.
  --gc_content_max_percent GC_CONTENT_MAX_PERCENT
                        Maximum GC percent for percent-mode clipping/axis.
  --gc_content_tick_interval GC_CONTENT_TICK_INTERVAL
                        GC content percent-mode large tick interval; alias for
                        --gc_content_large_tick_interval.
  --gc_content_large_tick_interval GC_CONTENT_LARGE_TICK_INTERVAL
                        GC content percent-mode large tick interval.
  --gc_content_small_tick_interval GC_CONTENT_SMALL_TICK_INTERVAL
                        GC content percent-mode small tick interval.
  --gc_content_tick_font_size GC_CONTENT_TICK_FONT_SIZE
                        GC content percent-mode tick label font size.
  --show_gc_content_axis, --hide_gc_content_axis
                        Show or hide the GC content percent-mode axis.
  --show_gc_content_ticks, --hide_gc_content_ticks
                        Show or hide GC content percent-mode ticks and labels.
  --show_skew           plot GC skew below genome (default: False).
  --align_center        Align genomes to the center (default: False).
  --evalue EVALUE       evalue threshold (default=1e-2)
  --bitscore BITSCORE   bitscore threshold (default=50)
  --identity IDENTITY   identity threshold (default=0)
  --alignment_length ALIGNMENT_LENGTH
                        minimum BLAST alignment length threshold (default=0)
  --pairwise_match_style {ribbon,curve}
                        Pairwise comparison link style: ribbon keeps straight
                        filled ribbons; curve draws curved filled ribbons that
                        preserve alignment spans.
  -k, --features FEATURES
                        Comma-separated list of feature keys to draw (default:
                        CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)
  --feature_shape TYPE=SHAPE
                        Feature shape override (repeatable): TYPE=SHAPE where
                        SHAPE is arrow or rectangle.
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "gray")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (optional; float; default: 2 pt for
                        genomes <= 50 kb, 0 pt for genomes >= 50 kb)
  --axis_stroke_color AXIS_STROKE_COLOR
                        Axis stroke color (str; default: auto: "lightgray", or
                        "dimgray" with --ruler_on_axis)
  --axis_stroke_width AXIS_STROKE_WIDTH
                        Axis stroke width (optional; float; default: 5 pt for
                        genomes <= 50 kb, 2 pt for genomes >= 50 kb)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (optional; str; default:
                        "lightgray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (optional; float; default: 5 pt for
                        genomes <= 50 kb, 1 pt for genomes >= 50 kb)
  --definition_font_size DEFINITION_FONT_SIZE
                        Definition font size (optional; float; default: 24 pt
                        for genomes <= 50 kb, 10 pt for genomes >= 50 kb)
  --plot_title PLOT_TITLE
                        Shared plot title text (optional).
  --plot_title_position {center,top,bottom}
                        Shared plot title position ("center", "top", "bottom";
                        default: "bottom").
  --plot_title_font_size PLOT_TITLE_FONT_SIZE
                        Shared plot title font size (optional; float; default:
                        32).
  --record_label RECORD_LABEL
                        Override record definition label (repeatable; order
                        matches input records)
  --label_font_size LABEL_FONT_SIZE
                        Label font size (optional; default: 24 pt for genomes
                        <= 50 kb, 5 pt for genomes >= 50 kb)
  --label_placement {auto,above_feature}
                        Linear label placement mode ("auto" or
                        "above_feature"; default: "auto"). "above_feature"
                        draws labels above features (or below negative-strand
                        features when --separate_strands is used).
  --label_rendering {auto,embedded_only,external_only}
                        Label rendering policy. "auto" embeds labels that fit
                        and routes the rest externally; "embedded_only" drops
                        external labels; "external_only" forces labels outside
                        feature bodies. Non-auto values cannot be combined with
                        --label_placement above_feature.
  --label_rotation LABEL_ROTATION
                        Linear label rotation in degrees (optional; float;
                        default: 0). In above_feature mode, rotated labels
                        start from the feature midpoint.
  --track_layout {above,middle,below}
                        Linear track layout mode ("above", "middle", or
                        "below"; default: "middle"). Aliases: "spreadout" ->
                        "above", "tuckin" -> "below".
  --track_axis_gap AUTO|PX
                        Gap between axis and nearest feature edge in pixels
                        for above/below layouts. Use 'auto' to derive it from
                        feature height.
  --linear_track_order LINEAR_TRACK_ORDER
                        Linear custom track shortcut order, for example
                        features,depth,gc_content,gc_skew.
  --linear_track_slot SLOT
                        Linear custom track slot:
                        <slot_id>:<renderer>@key=value,key=value. Repeat to
                        add slots.
  --linear_track_axis_index LINEAR_TRACK_AXIS_INDEX
                        Axis boundary index for linear custom track slots.
  --ruler_on_axis       Use each record axis as the ruler in linear mode.
                        Effective only with --scale_style ruler and
                        --track_layout above|below.
  -f, --format FORMAT   Comma-separated list of output file formats (svg,
                        interactive-svg, png, pdf, eps, ps; default: svg).
                        PNG/PDF/EPS/PS require CairoSVG.
  -l, --legend LEGEND   Legend position (default: "right"; "right", "left",
                        "top", "bottom", "none")
  --show_labels [{all,first,none}]
                        Show labels: no argument or 'all' (all records),
                        'first' (first record only), 'none' (no labels).
                        Default: 'none'
  --resolve_overlaps    Resolve overlaps (experimental; default: False).
  --label_whitelist LABEL_WHITELIST
                        Path to a TSV file for label whitelisting by regex
                        pattern (optional); mutually exclusive with
                        --label_blacklist
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label
                        blacklisting (optional); mutually exclusive with
                        --label_whitelist
  --qualifier_priority QUALIFIER_PRIORITY
                        Path to a TSV file defining qualifier priority for
                        labels (optional)
  --label_table LABEL_TABLE
                        Path to a TSV file defining post-filter label text
                        overrides (optional)
  --feature_table FEATURE_TABLE
                        Path to a TSV file defining per-feature visibility
                        overrides (optional)
  --feature_height FEATURE_HEIGHT
                        Feature vertical width (optional; float; default: 80
                        (pixels, 96 dpi) for genomes <= 50 kb, 20 for genomes
                        >= 50 kb)
  --gc_height GC_HEIGHT
                        GC content/skew vertical width (optional; float;
                        default: 20 (pixels, 96 dpi))
  --depth DEPTH [DEPTH ...]
                        Depth TSV file(s) in samtools depth format. Provide
                        one file for all records or one file per input record.
  --depth_track DEPTH [DEPTH ...]
                        Repeatable logical depth track. Each group accepts one
                        shared file or one file per input record.
  --depth_track_label LABEL [LABEL ...]
                        Depth track label(s). Provide one label or one per
                        --depth_track.
  --depth_track_color COLOR [COLOR ...]
                        Depth track fill color(s). Provide one color or one
                        per --depth_track.
  --show_depth          Show depth coverage track. Required only when no depth
                        file option is supplied.
  --depth_height DEPTH_HEIGHT
                        Depth track height for linear mode (in px; must be
                        > 0).
  --comparison_height COMPARISON_HEIGHT
                        Comparison block height (optional; float; optional;
                        default: 60 (pixels, 96 dpi))
  --scale_style {bar,ruler}
                        Style for the length scale (default: "bar"; "bar",
                        "ruler")
  --scale_stroke_color SCALE_STROKE_COLOR
                        Scale bar/ruler stroke color (optional; str; default:
                        "black")
  --scale_stroke_width SCALE_STROKE_WIDTH
                        Scale bar/ruler stroke width (optional; float;
                        default: 3 (pt))
  --scale_font_size SCALE_FONT_SIZE
                        Scale bar/ruler font size (optional; float; default:
                        24 (pt) for genomes <= 50 kb, 16 for genomes >= 50
                        kb).
  --ruler_label_font_size RULER_LABEL_FONT_SIZE
                        Ruler label font size (optional; float). Overrides
                        --scale_font_size when both are set.
  --ruler_label_color RULER_LABEL_COLOR
                        Ruler label color (optional; str; default follows axis
                        color when --ruler_on_axis is active, otherwise
                        black).
  --scale_interval SCALE_INTERVAL
                        Manual tick interval for "ruler" scale style (in bp).
                        Overrides automatic calculation; optional
  --legend_box_size LEGEND_BOX_SIZE
                        Legend box size (optional; float; default: 24 (pixels,
                        96 dpi) for genomes <= 50 kb, 20 for genomes >= 50
                        kb).
  --legend_font_size LEGEND_FONT_SIZE
                        Legend font size (optional; float; default: 20 (pt)
                        for genomes <= 50 kb, 16 for genomes >= 50 kb).
  --normalize_length    Normalize record length (experimental; default:
                        False).
  --region REGION       Crop a region (repeatable). Format: record_id:start-
                        end[:rc], #index:start-end[:rc], or
                        file:record_selector:start-end[:rc]. Coordinates are
                        1-based inclusive. For multiple records without
                        selectors, provide one spec per record in input order
                        (file order, then record order within each file).
  --record_id RECORD_ID
                        Select a record by ID or #index per input file
                        (repeatable; order matches input files). Use an empty
                        value to skip selection for a file.
  --reverse_complement REVERSE_COMPLEMENT
                        Reverse complement record per input file (repeatable;
                        order matches input files). Accepted values: 1/0,
                        true/false, yes/no.
```

For `--protein_blastp_mode`, gbdraw first uses a bundled native LOSAT binary
when one is available. The current package bundles LOSAT for Linux x86_64.
macOS and Windows packages do not currently include bundled LOSAT binaries; if
no native LOSAT executable is available, install NCBI BLAST+ and make `blastp`
available on `PATH`, or pass it explicitly with `--ncbi_blastp_bin`. You can
still force a native LOSAT executable on any platform with `--losatp_bin`.
NCBI BLAST+ fallback produces compatible outfmt 6 protein comparisons, but its
hit set is not guaranteed to be identical to LOSAT.

## Related Guides

- [Quickstart](./QUICKSTART.md)
- [Recipes](./RECIPES.md)
- [Tutorials](./TUTORIALS/TUTORIALS.md)

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | **CLI Reference** | [Gallery](./GALLERY.md) | [FAQ](./FAQ.md) | [About](./ABOUT.md)
