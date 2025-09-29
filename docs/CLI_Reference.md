
# Usage

```bash
$ gbdraw -h
gbdraw v. 0.5.0: A diagram generator for small genomes

Usage:
  gbdraw <subcommand> [options]

Subcommands:
  circular  Generate a circular genome diagram
  linear    Generate a linear genome diagram
  gui       Launch the graphical user interface (requires Streamlit)

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

Additional Information:
  - For full documentation, visit: https://github.com/satoshikawato/gbdraw/
  - For issues and source code, visit the GitHub repository: https://github.com/satoshikawato/gbdraw/
  - For support, contact: kawato[at]kaiyodai.ac.jp
```
## Circular mode
```bash
$ gbdraw circular -h
usage: gbdraw [-h] [--gbk [GBK_FILE ...]] [--gff [GFF3_FILE ...]] [--fasta [FASTA_FILE ...]] [-o OUTPUT] [-p PALETTE]
              [-t TABLE] [-d DEFAULT_COLORS] [-n NT] [-w WINDOW] [-s STEP] [--species SPECIES] [--strain STRAIN]
              [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR] [--block_stroke_width BLOCK_STROKE_WIDTH]
              [--axis_stroke_color AXIS_STROKE_COLOR] [--axis_stroke_width AXIS_STROKE_WIDTH]
              [--line_stroke_color LINE_STROKE_COLOR] [--line_stroke_width LINE_STROKE_WIDTH]
              [--definition_font_size DEFINITION_FONT_SIZE] [--label_font_size LABEL_FONT_SIZE] [-f FORMAT]
              [--suppress_gc] [--suppress_skew] [-l LEGEND] [--separate_strands] [--track_type TRACK_TYPE]
              [--show_labels] [--allow_inner_labels] [--label_whitelist LABEL_WHITELIST]
              [--label_blacklist LABEL_BLACKLIST] [--qualifier_priority QUALIFIER_PRIORITY]
              [--outer_label_x_radius_offset OUTER_LABEL_X_RADIUS_OFFSET]
              [--outer_label_y_radius_offset OUTER_LABEL_Y_RADIUS_OFFSET]
              [--inner_label_x_radius_offset INNER_LABEL_X_RADIUS_OFFSET]
              [--inner_label_y_radius_offset INNER_LABEL_Y_RADIUS_OFFSET]

Generate genome diagrams in PNG/PDF/SVG/PS/EPS. Diagrams for multiple entries are saved separately.

options:
  -h, --help            show this help message and exit
  --gbk [GBK_FILE ...]  Genbank/DDBJ flatfile
  --gff [GFF3_FILE ...]
                        GFF3 file (instead of --gbk; --fasta is required)
  --fasta [FASTA_FILE ...]
                        FASTA file (required with --gff)
  -o, --output OUTPUT   output file prefix (default: accession number of the sequence)
  -p, --palette PALETTE
                        Palette name (default: default)
  -t, --table TABLE     color table (optional)
  -d, --default_colors DEFAULT_COLORS
                        TSV file that overrides the color palette (optional)
  -n, --nt NT           dinucleotide (default: GC).
  -w, --window WINDOW   window size (default: 1000)
  -s, --step STEP       step size (default: 100)
  --species SPECIES     Species name (optional; e.g. "<i>Escherichia coli</i>", "<i>Ca.</i> Hepatoplasma
                        crinochetorum")
  --strain STRAIN       Strain/isolate name (optional; e.g. "K-12", "Av")
  -k, --features FEATURES
                        Comma-separated list of feature keys to draw (default:
                        CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "gray")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (float; default: 0)
  --axis_stroke_color AXIS_STROKE_COLOR
                        Axis stroke color (str; default: "gray")
  --axis_stroke_width AXIS_STROKE_WIDTH
                        Axis stroke width (float; default: 1.0)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (str; default: "gray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (float; default: 1.0)
  --definition_font_size DEFINITION_FONT_SIZE
                        Definition font size (optional; default: 18)
  --label_font_size LABEL_FONT_SIZE
                        Label font size (optional; default: 16 for short genomes, 8 for long genomes)
  -f, --format FORMAT   Comma-separated list of output file formats (default: png)
  --suppress_gc         Suppress GC content track (default: False).
  --suppress_skew       Suppress GC skew track (default: False).
  -l, --legend LEGEND   Legend position (default: "right"; "left", "right", "upper_left", "upper_right", "lower_left",
                        "lower_right", "none")
  --separate_strands    Separate strands (default: False).
  --track_type TRACK_TYPE
                        Track type (default: "tuckin"; "tuckin", "middle", "spreadout")
  --show_labels         Show feature labels (default: False).
  --allow_inner_labels  Place labels inside the circle (default: False). If enabled, labels are placed both inside and
                        outside the circle, and gc and skew tracks are not shown.
  --label_whitelist LABEL_WHITELIST
                        path to a file for label whitelisting (optional); mutually exclusive with --label_blacklist
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label blacklisting (optional); mutually
                        exclusive with --label_whitelist
  --qualifier_priority QUALIFIER_PRIORITY
                        Path to a TSV file defining qualifier priority for labels (optional)
  --outer_label_x_radius_offset OUTER_LABEL_X_RADIUS_OFFSET
                        Outer label x-radius offset factor (float; default from config)
  --outer_label_y_radius_offset OUTER_LABEL_Y_RADIUS_OFFSET
                        Outer label y-radius offset factor (float; default from config)
  --inner_label_x_radius_offset INNER_LABEL_X_RADIUS_OFFSET
                        Inner label x-radius offset factor (float; default from config)
  --inner_label_y_radius_offset INNER_LABEL_Y_RADIUS_OFFSET
                        Inner label y-radius offset factor (float; default from config)
```

## Linear mode

```bash
$ gbdraw linear -h
usage: gbdraw [-h] [--gbk [GBK_FILE ...]] [--gff [GFF3_FILE ...]] [--fasta [FASTA_FILE ...]] [-b [BLAST ...]]
              [-t TABLE] [-p PALETTE] [-d DEFAULT_COLORS] [-o OUTPUT] [-n NT] [-w WINDOW] [-s STEP]
              [--separate_strands] [--show_gc] [--align_center] [--evalue EVALUE] [--bitscore BITSCORE]
              [--identity IDENTITY] [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR]
              [--block_stroke_width BLOCK_STROKE_WIDTH] [--axis_stroke_color AXIS_STROKE_COLOR]
              [--axis_stroke_width AXIS_STROKE_WIDTH] [--line_stroke_color LINE_STROKE_COLOR]
              [--line_stroke_width LINE_STROKE_WIDTH] [--definition_font_size DEFINITION_FONT_SIZE]
              [--label_font_size LABEL_FONT_SIZE] [-f FORMAT] [-l LEGEND] [--show_labels] [--resolve_overlaps]
              [--label_whitelist LABEL_WHITELIST | --label_blacklist LABEL_BLACKLIST]
              [--qualifier_priority QUALIFIER_PRIORITY]

Generate plot in PNG/PDF/SVG/PS/EPS.

options:
  -h, --help            show this help message and exit
  --gbk [GBK_FILE ...]  Genbank/DDBJ flatfile
  --gff [GFF3_FILE ...]
                        GFF3 file (instead of --gbk; --fasta is required)
  --fasta [FASTA_FILE ...]
                        FASTA file (required with --gff)
  -b, --blast [BLAST ...]
                        input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)
  -t, --table TABLE     color table (optional)
  -p, --palette PALETTE
                        Palette name (default: default)
  -d, --default_colors DEFAULT_COLORS
                        TSV file that overrides the color palette (optional)
  -o, --output OUTPUT   output file prefix (default: out)
  -n, --nt NT           dinucleotide skew (default: GC).
  -w, --window WINDOW   window size (default: 1000)
  -s, --step STEP       step size (default: 100)
  --separate_strands    separate forward and reverse strands (default: False). Features of undefined strands are shown
                        on the forward strand.
  --show_gc             plot GC content below genome (default: False).
  --align_center        Align genomes to the center (default: False).
  --evalue EVALUE       evalue threshold (default=1e-2)
  --bitscore BITSCORE   bitscore threshold (default=50)
  --identity IDENTITY   identity threshold (default=0)
  -k, --features FEATURES
                        Comma-separated list of feature keys to draw (default:
                        CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,repeat_region)
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "black")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (float; default: 0)
  --axis_stroke_color AXIS_STROKE_COLOR
                        Axis stroke color (str; default: "lightgray")
  --axis_stroke_width AXIS_STROKE_WIDTH
                        Axis stroke width (float; default: 2.0)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (str; default: "gray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (float; default: 1.0)
  --definition_font_size DEFINITION_FONT_SIZE
                        Definition font size (optional; default: 10)
  --label_font_size LABEL_FONT_SIZE
                        Label font size (optional; default: 16 for short genomes, 5 for long genomes)
  -f, --format FORMAT   Comma-separated list of output file formats (default: png)
  -l, --legend LEGEND   Legend position (default: "right"; "right", "left", "none")
  --show_labels         Show labels
  --resolve_overlaps    Resolve overlaps (experimental; default: False).
  --label_whitelist LABEL_WHITELIST
                        path to a file for label whitelisting (optional); mutually exclusive with --label_blacklist
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label blacklisting (optional); mutually
                        exclusive with --label_whitelist
  --qualifier_priority QUALIFIER_PRIORITY
                        Path to a TSV file defining qualifier priority for labels (optional)
```
