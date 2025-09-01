[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/version.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/platforms.svg)](https://anaconda.org/bioconda/gbdraw)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/gbdraw/README.html)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/latest_release_date.svg)](https://anaconda.org/bioconda/gbdraw)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/gbdraw/badges/license.svg)](https://anaconda.org/bioconda/gbdraw)
[![Ask DeepWiki](https://deepwiki.com/badge.svg)](https://deepwiki.com/satoshikawato/gbdraw)

# gbdraw
![gbdraw](https://github.com/satoshikawato/gbdraw/blob/main/examples/gbdraw_social_preview.png)
`gbdraw` is a command-line/GUI tool designed for creating detailed diagrams of microbial genomes. 
`gbdraw` accepts GenBank/EMBL/DDBJ-format annotated genomes as input and outputs a visual representation of the genomes in SVG/PNG/PDF/EPS/PS formats.

**Try gbdraw Web App!** [https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

**Try gbdraw on Colab Notebook!** [![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)

## Features
- Circular and linear diagrams: Generates both circular and linear representations of genome structures.
- Customizable inputs: Supports Genbank/DDBJ flat files with options for color customization.
- Various output formats: Vector and raster graphics suitable for publication and further editing.

## Dependencies
- [Python](https://www.python.org/) >=3.10
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [svgwrite](https://github.com/mozman/svgwrite)
- [CairoSVG](https://cairosvg.org/)
- [Liberation Fonts](https://github.com/liberationfonts/liberation-fonts) (bundled; SIL Open Font_License 1.1)

## Use without local installation
### Streamlit Web App
A GUI web app of `gbdraw` (latest commit on `main` branch) is available on Streamlit without any local installation:
[https://gbdraw.streamlit.app/](https://gbdraw.streamlit.app/)

### Colab Notebook (Google account required)
You can try `gbdraw` (latest release) on Google Colaboratory without any local installation:

[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/satoshikawato/gbdraw/blob/main/gbdraw_colab.ipynb)

## Local Installation
**Prerequisite:** Make sure you have a [conda](https://docs.conda.io/en/latest/)-compatible package manager—[mamba](https://github.com/mamba-org/mamba) ,[micromamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html), [miniforge](https://github.com/conda-forge/miniforge) or plain conda—already installed and on your `$PATH`. All steps below assume you run the commands in such an environment.
### Bioconda (recommended)
`gbdraw` is available on the Bioconda channel.
```bash
conda create -n gbdraw-0.4.0 -y -c conda-forge -c bioconda gbdraw=0.4.0
conda activate gbdraw-0.4.0
```
### Local build (development version)
To use the latest development version, clone the repository yourself using `git` and build the package locally with [conda-build](https://anaconda.org/anaconda/conda-build).
```bash
# 1. Clone the source
git clone https://github.com/satoshikawato/gbdraw.git
cd gbdraw/

# 2. Make sure conda-build is installed
mamba install -y conda-build             # or: conda install conda-build

# 3. Build the package locally
conda-build .

# 4. Create an isolated environment from the locally built package
mamba create -n gbdraw -y  -c conda-forge -c bioconda -c local gbdraw

# 5. Activate the environment
mamba activate gbdraw
```


## Usage
```bash
$ gbdraw 
gbdraw v. 0.4.0: A diagram generator for small genomes

Usage:
  gbdraw <subcommand> [options]

Subcommands:
  circular  Generate a circular genome diagram
  linear    Generate a linear genome diagram

For each subcommand, you can get additional help by running:
  gbdraw <subcommand> --help

Examples:
  gbdraw circular -i input.gb
  gbdraw linear -i input.gb

Options (examples):
  -i, --input          Input GenBank file(s) (required)
  -o, --output         Output file prefix (optional)
  -t, --table          Color table file (optional)
  -b, --blast          BLAST result file in tab-separated format (-outfmt 6 or 7) (optional; currently implemented for linear mode only)

Additional Information:
  - For full documentation, visit: https://github.com/satoshikawato/gbdraw/
  - For issues and source code, visit the GitHub repository: https://github.com/satoshikawato/gbdraw/
  - For support, contact: kawato[at]kaiyodai.ac.jp
```
### Circular genome
```bash
$ gbdraw circular -h
usage: gbdraw [-h] -i [INPUT ...] [-o OUTPUT] [-p PALETTE] [-t TABLE] [-d DEFAULT_COLORS] [-n NT] [-w WINDOW]
              [-s STEP] [--species SPECIES] [--strain STRAIN] [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR]
              [--block_stroke_width BLOCK_STROKE_WIDTH] [--axis_stroke_color AXIS_STROKE_COLOR]
              [--axis_stroke_width AXIS_STROKE_WIDTH] [--line_stroke_color LINE_STROKE_COLOR]
              [--line_stroke_width LINE_STROKE_WIDTH] [--definition_font_size DEFINITION_FONT_SIZE]
              [--label_font_size LABEL_FONT_SIZE] [-f FORMAT] [--suppress_gc] [--suppress_skew] [-l LEGEND]
              [--separate_strands] [--track_type TRACK_TYPE] [--show_labels] [--allow_inner_labels]
              [--label_blacklist LABEL_BLACKLIST] [--qualifier_priority QUALIFIER_PRIORITY]
              [--outer_label_x_radius_offset OUTER_LABEL_X_RADIUS_OFFSET]
              [--outer_label_y_radius_offset OUTER_LABEL_Y_RADIUS_OFFSET]
              [--inner_label_x_radius_offset INNER_LABEL_X_RADIUS_OFFSET]
              [--inner_label_y_radius_offset INNER_LABEL_Y_RADIUS_OFFSET]

Generate genome diagrams in PNG/PDF/SVG/PS/EPS. Diagrams for multiple entries are saved separately.

options:
  -h, --help            show this help message and exit
  -i, --input [INPUT ...]
                        Genbank/DDBJ flatfile (required)
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
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label blacklisting (optional)
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
#### <i>Haemophilus influenzae</i>
`gbdraw` automatically identifies and displays the organism and strain name from the sequence record. However, these names are not italicized by default. For example:
```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/931/575/GCF_000931575.1_ASM93157v1/GCF_000931575.1_ASM93157v1_genomic.gbff.gz # download genome
gunzip GCF_000931575.1_ASM93157v1_genomic.gbff.gz # decompress .gz file
gbdraw circular -i GCF_000931575.1_ASM93157v1_genomic.gbff -o Haemophilus_influenzae -f svg
```
![hinfluenzae](https://github.com/satoshikawato/gbdraw/blob/main/examples/Haemophilus_influenzae.svg)
#### <i>Escherichia coli</i> K-12
To italicize a portion of the organism name, you can use the <i></i> tags in the --species and --strain parameters. This will format the specified text in italics. The following command will render the species name "_Escherichia coli_" in italics, while keeping "K-12" in standard text (the organim name will be overridden):
```bash
gbdraw circular -i NC_000913.gb --species "<i>Escherichia coli</i>" --strain "K-12" -f svg --separate_strands
```
![ecoli](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000913.svg)

#### <i>Vibrio cholerae</i> Strain O395 (GCF_000016245.1)
For GenBank files containing multiple entries, `gbdraw` saves each entry as a separate file. Here's how to do it for the _Vibrio cholerae_ strain O395 genome, which has two chromosomes:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/245/GCF_000016245.1_ASM1624v1/GCF_000016245.1_ASM1624v1_genomic.gbff.gz # download genome
gunzip GCF_000016245.1_ASM1624v1_genomic.gbff.gz # decomperss .gz file
gbdraw circular -i GCF_000016245.1_ASM1624v1_genomic.gbff --species "<i>Vibrio cholerae</i>" --strain "O395" -f svg --separate_strands --track_type middle # Draw genome; results in "NC_009457.svg" for Chromosome I and "NC_009456.svg" for Chromosome II

```
<i>Vibrio cholerae</i> Chromosome I
![Vibrio cholerae chromosome I](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009457.svg)
<i>Vibrio cholerae</i> Chromosome II
![Vibrio cholerae chromosome II](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009456.svg)
### Linear genome
`gbdraw linear`
```bash
$ gbdraw linear -h
usage: gbdraw [-h] -i [INPUT ...] [-b [BLAST ...]] [-t TABLE] [-p PALETTE] [-d DEFAULT_COLORS] [-o OUTPUT] [-n NT]
              [-w WINDOW] [-s STEP] [--separate_strands] [--show_gc] [--align_center] [--evalue EVALUE]
              [--bitscore BITSCORE] [--identity IDENTITY] [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR]
              [--block_stroke_width BLOCK_STROKE_WIDTH] [--axis_stroke_color AXIS_STROKE_COLOR]
              [--axis_stroke_width AXIS_STROKE_WIDTH] [--line_stroke_color LINE_STROKE_COLOR]
              [--line_stroke_width LINE_STROKE_WIDTH] [--definition_font_size DEFINITION_FONT_SIZE]
              [--label_font_size LABEL_FONT_SIZE] [-f FORMAT] [-l LEGEND] [--show_labels] [--resolve_overlaps]
              [--label_blacklist LABEL_BLACKLIST] [--qualifier_priority QUALIFIER_PRIORITY]

Generate plot in PNG/PDF/SVG/PS/EPS.

options:
  -h, --help            show this help message and exit
  -i, --input [INPUT ...]
                        genbank (required)
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
  --label_blacklist LABEL_BLACKLIST
                        Comma-separated keywords or path to a file for label blacklisting (optional)
  --qualifier_priority QUALIFIER_PRIORITY
                        Path to a TSV file defining qualifier priority for labels (optional)
```
### Human herpesvirus 6 (HHV-6)
`gbdraw` can draw linear genomes and pairwise matches depicting similar genomic regions.
Perform a pariwise BLASTN analysis between the genomes of HHV-6A and HHV-6B:
```bash
blastn -query NC_000898.fasta -subject NC_001664.fasta -outfmt 7 -out NC_000898_NC_001664.blastn.out
```
Use `gbdraw linear` to visualize the comparison, aligning genomes to the center and separating forward and reverse strands:
```bash
umamba activate blast-2.16.0
blastn -query NC_000898.fasta -subject NC_001664.fasta -outfmt 7 -out NC_000898_NC_001664.blastn.out
umamba deactivate
umamba activate gbdraw-0.1.0
gbdraw linear -i NC_000898.gb NC_001664.gb -b NC_000898_NC_001664.blastn.out --resolve_overlaps --align_center --separate_strands -o HHV-6 -f svg --block_stroke_width 0.5

umamba deactivate
```
![HHV-6](https://github.com/satoshikawato/gbdraw/blob/main/examples/HHV-6.svg)

### <i>Candidatus</i> Hepatoplasmataceae (mollicutes)
```bash
tblastx -query Fukuoka2020.fasta -subject Av-JP.fasta  -outfmt 7 -out Fukuoka2020_Av-JP.tblastx.out
tblastx -query Av-JP.fasta -subject Ps-JP.fasta  -outfmt 7 -out Av-JP_Ps-JP.tblastx.out
tblastx -query Ps-JP.fasta -subject Tokyo2021.fasta  -outfmt 7 -out Ps-JP_Tokyo2021.tblastx.out
tblastx -query Tokyo2021.fasta -subject Av.fasta  -outfmt 7 -out Tokyo2021_Av.tblastx.out
gbdraw linear -i Fukuoka2020.gb Av-JP.gb Ps-JP.gb Tokyo2021.gb Av.gb -b Fukuoka2020_Av-JP.tblastx.out  Av-JP_Ps-JP.tblastx.out  Ps-JP_Tokyo2021.tblastx.out  Tokyo2021_Av.tblastx.out -o hepatoplasmataceae --align_center --bitscore 50 --evalue 1e-3 --separate_strands -f svg
```
![hepatoplasmataceae](https://github.com/satoshikawato/gbdraw/blob/main/examples/hepatoplasmataceae.svg)
## Advanced customization

### Color palettes

`gbdraw` ships with [a total of 55 color palettes](https://github.com/satoshikawato/gbdraw/blob/main/examples/color_palette_examples.md). Choose a palette with **`-p/--palette <name>`** or override individual colours via TSV files.
![palettes_combined_image_1.png](https://github.com/satoshikawato/gbdraw/blob/main/examples/palettes_combined_image_1.png)

#### Examples
##### autumn
![autumn](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_autumn.svg)
##### forest
![forest](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_forest.svg)
##### fugaku
![fugaku](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_fugaku.svg)
##### lavender_fields
![lavender_fields](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_lavender_fields.svg)
##### matcha_whispers
![matcha_whispers](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_matcha_whispers.svg)
##### sakura
![sakura](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_sakura.svg)
##### tropical
![tropical](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_tropical.svg)
##### zen
![zen](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_tuckin_separate_strands_zen.svg)
See [this page](https://github.com/satoshikawato/gbdraw/blob/main/examples/color_palette_examples.md) for the examples of all 55 palettes.

### Customizing colors with configuration files
`gbdraw` supports two complementary mechanisms for overriding the default colours:

| Method | Purpose |
| ------ | ------- |
| **Default-override table** (`-d`) | Replace colours for *entire* feature classes |
| **Feature-specific table** (`-t`) | Colour *individual* features that match user-defined rules |

Both tables are tab separated. Colours may be given as any of the [147 color names defined by the SVG specification](https://johndecember.com/html/spec/colorsvg.html) or in [hexadecimal format](https://htmlcolorcodes.com/) (`#RRGGBB`).

#### Overriding the default color table
The following `modified_default_colors.tsv` turns CDS gray (`#d3d3d3`). Other features remain the same as default:
```modified_default_colors.tsv
CDS	#d3d3d3
```

| feature type | color |
| ------ | ------- |
| CDS | #d3d3d3 |

#### Feature-specific color table
For more precise control, especially for CDS features, you can specify colors based on specific attributes like `product` or `note`.
This `custom_color_table.tsv` color only those CDS whose `product` qualifier matches the given regex:
```custom_color_table.tsv
CDS	product	wsv.*-like protein	#47b8f8	WSSV-like proteins
CDS	product	baculoviral IAP repeat-containing protein	yellow	BIRP
CDS	product	tyrosine recombinase	red	tyrosine recombinase
```
This means:
| feature type | target qualifier | qualifier value regex (Python) | color | legend text |
| ------ | ------- | ------- | ------- | ------- |
| CDS | product | wsv.*-like protein | #47b8f8 | WSSV-like proteins |
| CDS | product | baculoviral IAP repeat-containing protein | yellow | BIRP |
| CDS | product | tyrosine recombinase | red | tyrosine recombinase |

```bash
gbdraw circular -i LC738868.gb -o LC738868_middle_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv 
```
![MjeNMV](https://github.com/satoshikawato/gbdraw/blob/main/examples/LC738868_middle_separate_strands.png)

### Multiple genome alignments
`gbdraw linear` can also be used for visualizing multi-genome alignments, providing a comparative view of different genomes. This feature is particularly useful for identifying conserved regions and variations across multiple genomes.
**NOTE:** Currently the BLAST output files must be specified in exactly the same order as the genbank files.
Example command for a linear genome diagram with multi-genome alignment:
```bash
# Pairwise TBLASTX search 
tblastx -query MjeNMV.fasta -subject MelaMJNV.fasta -outfmt 7 -out MjeNMV.MelaMJNV.tblastx.out 
tblastx -query MelaMJNV.fasta -subject PemoMJNVA.fasta -outfmt 7 -out MelaMJNV.PemoMJNVA.tblastx.out 
tblastx -query PemoMJNVA.fasta -subject PeseMJNV.fasta -outfmt 7 -out PemoMJNVA.PeseMJNV.tblastx.out 
tblastx -query PeseMJNV.fasta -subject PemoMJNVB.fasta -outfmt 7 -out PeseMJNV.PemoMJNVB.tblastx.out 
tblastx -query PemoMJNVB.fasta -subject LvMJNV.fasta -outfmt 7 -out PemoMJNVB.LvMJNV.tblastx.out 
tblastx -query LvMJNV.fasta -subject TrcuMJNV.fasta -outfmt 7 -out LvMJNV.TrcuMJNV.tblastx.out 
tblastx -query TrcuMJNV.fasta -subject MellatMJNV.fasta -outfmt 7 -out TrcuMJNV.MellatMJNV.tblastx.out 
tblastx -query MellatMJNV.fasta -subject MeenMJNV.fasta -outfmt 7 -out MellatMJNV.MeenMJNV.tblastx.out 
tblastx -query MeenMJNV.fasta -subject MejoMJNV.fasta -outfmt 7 -out MeenMJNV.MejoMJNV.tblastx.out 

# gbdraw
gbdraw linear \
-i \
./in_gbk/MjeNMV.gb \
./in_gbk/MelaMJNV.gb \
./in_gbk/PemoMJNVA.gb \
./in_gbk/PeseMJNV.gb \
./in_gbk/PemoMJNVB.gb \
./in_gbk/LvMJNV.gb \
./in_gbk/TrcuMJNV.gb \
./in_gbk/MetlamMJNV.gb \
./in_gbk/MeenMJNV.gb \
./in_gbk/MejoMJNV.gb \
-b \
./in_fna/MjeNMV.MelaMJNV.tblastx.out \
./in_fna/MelaMJNV.PemoMJNVA.tblastx.out \
./in_fna/PemoMJNVA.PeseMJNV.tblastx.out \
./in_fna/PeseMJNV.PemoMJNVB.tblastx.out \
./in_fna/PemoMJNVB.LvMJNV.tblastx.out \
./in_fna/LvMJNV.TrcuMJNV.tblastx.out \
./in_fna/TrcuMJNV.MetlamMJNV.tblastx.out \
./in_fna/MetlamMJNV.MeenMJNV.tblastx.out \
./in_fna/MeenMJNV.MejoMJNV.tblastx.out \
-t majani_custom_color_table.tsv \
-d modified_default_colors.tsv \
--block_stroke_width 1 \
--block_stroke_color gray \
--align_center \
--separate_strands \
-o majani -f svg
```
![majaniviruses](https://github.com/satoshikawato/gbdraw/blob/main/examples/majani.svg)

### Feature labels
`gbdraw` can draw feature labels. By default, most genic features have the value of `product` qualifer as the label text.
| feature type | priority of qualifiers used for the label |
| ------ | ------- |
| 'CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA', 'misc_RNA', 'gene' | 'product', 'gene', 'note' |
| repeat_region |'rpt_family', 'note'|
| other features | note |
```bash
gbdraw circular -i AP027280.gb -f svg --block_stroke_width 1 --block_stroke_color gray --track_type spreadout --show_labels
```
![WSSV](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027280.svg)
```bash
gbdraw circular -i NC_012920.gb -f svg --block_stroke_width 2 --block_stroke_color gray  --show_labels -w 100 -s 10
```
![HsmtDMA](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_012920.svg)

#### Label size adjustment
The font size of the Label text can be adjusted via `--label_font_size` option.
#### Label text blacklist 
You can selectively hide some labels via `--label_blacklist` option. For example, supressing `hypothetical protein` and other less-informative labels reduces label redundancy and makes it easier to visualize larger genomes without layout issues.
`--label_blacklist` accepts comma-separated keywords or path to a file containing single keyword per line.
```bash
gbdraw circular -i AP027280.gb \
-f svg --block_stroke_width 1 \
--track_type middle \
--show_labels \
--label_blacklist hypothetical \
--outer_label_x_radius_offset 0.95 \
--outer_label_y_radius_offset 0.95 \
-o AP027280_middle_blacklist_ox0.95_oy0.95
```
![AP027280_middle_blacklist_ox0.95_oy0.95](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027280_middle_blacklist_ox0.95_oy0.95.svg)

```bash
gbdraw circular \
-i AP027078.gb \
-f svg \
--track_type middle \
--separate_strands \
--show_labels \
--label_font_size 4 \
--label_blacklist hypothetical \
--outer_label_y_radius_offset 1.05 \
-o AP027078_middle_separate_strands_show_labels_label_blacklist_oy1.05
```
![AP027078_middle_separate_strands_show_labels_label_blacklist_oy1.05](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027078_middle_separate_strands_show_labels_label_blacklist_oy1.05.svg)

#### Feature qualifier priority
The following `qualifier_priority.tsv` designated by `--qualifier_priority` opton specifies which qualifier should be used for the label text of a given feature type. Other features remain the same as default:

```modified_default_colors.tsv
CDS gene
```
| feature type | qualifier |
| ------ | ------- |
| CDS | gene |

```bash
gbdraw circular \
-i NC_012920.gb \
-f svg --track_type middle \
--species "<i>Homo sapiens</i>" \
--block_stroke_width 2 \
--axis_stroke_width 5 \
--allow_inner_labels \
--show_labels \
--qualifier_priority qualifier_priority.tsv \
-o NC_012920_middle_qualifier_priority_inner_axis5_def28_italic \
--definition_font_size 28
```
![NC_012920_middle_qualifier_priority_inner_axis5_def28_italic](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_012920_middle_qualifier_priority_inner_axis5_def28_italic.svg)

## Planned features
- Feature label whitelist (display labels only for specified features)
- GFF3 support (will require additonal dependencies; GC content/skew will not be displayed unless the sequence is included in the GFF3 file or is provided separately as a FASTA file; species/strain must be explicitly designated by the user)
- PAF/MAF support (pairwise matches) 
- Caching (particularly relevant for the webapp; unsure how to impelement)
- Feature overlap resolution (overlapping genes, transcript isoforms etc.;unsure how to impelment mutiple tracks with label overlap resolution)
- Custom tracks (read depth, motifs, etc.)
## Known issues
- **Trans-introns** are not currently visualized.
- **Mixed-format text** (e.g., combining italic and block elements like `<i>Ca.</i> Tyloplasma litorale`) cannot be reliably converted from SVG to PDF/PNG/EPS/PS.  
  → As a workaround, export to **SVG format** and convert to other formats using external tools like [**Inkscape**](https://inkscape.org/).

## Bug reports and suggestions
Please feel free to submit a new issue if you find a bug or have a suggestion:
https://github.com/satoshikawato/gbdraw/issues
## Citation
If you found `gbdraw` useful, I would be very happy if you could cite the following URL:
https://github.com/satoshikawato/gbdraw/


## Note
This package was inspired by a variety of existing programs and tools in the bioinformatics field. The development of gbdraw was guided by the features and capabilities of these tools, and we acknowledge their influence in shaping our approach to genome visualization and analysis:

[CGView](https://cgview.ca/)

[ACT](https://www.sanger.ac.uk/tool/artemis-comparison-tool-act/)

[SnapGene Viewer](https://www.snapgene.com/snapgene-viewer)

[DNA Features Viewer](https://github.com/Edinburgh-Genome-Foundry/DnaFeaturesViewer)

[pyGenomeViz](https://github.com/moshi4/pyGenomeViz)

[Circos](https://www.circos.ca/)

[GenomeMatcher](https://www.ige.tohoku.ac.jp/joho/gmProject/gmhomeJP.html)

[GenomeDiagram](https://biopython-tutorial.readthedocs.io/en/latest/notebooks/17%20-%20Graphics%20including%20GenomeDiagram.html)

[GenoVi](https://github.com/robotoD/GenoVi)

[OrganellarGenomeDRAW (OGDRAW)](https://chlorobox.mpimp-golm.mpg.de/OGDraw.html)

The core functionality of gbdraw has evolved from [a set of Python scripts](https://github.com/satoshikawato/bio_small_scripts/) written back in 2022. Key among these are:

[plot_circular_genome.py](https://github.com/satoshikawato/bio_small_scripts/blob/main/plot_circular_genome.py)

[plot_linear_genome.py](https://github.com/satoshikawato/bio_small_scripts/blob/main/plot_linear_genome.py)




