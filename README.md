# gbdraw
![gbdraw](https://github.com/satoshikawato/gbdraw/blob/main/examples/gbdraw_preview.png)
`gbdraw` is a command-line tool designed for creating detailed diagrams of microbial genomes. 
`gbdraw` accepts GenBank/EMBL/DDBJ-format annotated genomes as input and outputs a visual representation of the genomes in SVG/PNG/PDF/EPS/PS formats.

**NOTE:** `gbdraw` is currently a work in progress and is actively under development, with the goal of releasing it as a user-friendly conda package in the near future. The current repository is a preview and does not yet host a fully functional release version. Stay tuned for updates and releases!
## Features
- Circular and linear diagrams: Generates both circular and linear representations of genome structures.
- Customizable inputs: Supports Genbank/DDBJ flat files with options for color customization.
- Various output formats: Vector and raster graphics suitable for publication and further editing.
## Dependencies
- [Python](https://www.python.org/) >=3.12
- [Biopython](https://biopython.org/)
- [pandas](https://pandas.pydata.org/)
- [svgwrite](https://github.com/mozman/svgwrite)
- [CairoSVG](https://cairosvg.org/)
- [Liberation Fonts](https://github.com/liberationfonts/liberation-fonts) (bundled; SIL Open Font_License 1.1)
## Usage
```bash
gbdraw -h
gbdraw v. 0.1.0: A diagram generator for small genomes

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
```
### Circular genome
```bash
gbdraw circular -h
usage: gbdraw [-h] -i [INPUT ...] [-o OUTPUT] [-t TABLE] [-d DEFAULT_COLORS] [-n NT] [-w WINDOW] [-s STEP] [--species SPECIES]
              [--strain STRAIN] [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR] [--block_stroke_width BLOCK_STROKE_WIDTH]
              [--line_stroke_color LINE_STROKE_COLOR] [--line_stroke_width LINE_STROKE_WIDTH] [-f FORMAT] [--suppress_gc]
              [--suppress_skew]

Generate genome diagrams in PNG/PDF/SVG/PS/EPS. Diagrams for multiple entries are saved separately (hence the lack of output file
name option).

options:
  -h, --help            show this help message and exit
  -i [INPUT ...], --input [INPUT ...]
                        Genbank/DDBJ flatfile (required)
  -o OUTPUT, --output OUTPUT
                        output file prefix (default: accession number of the sequence)
  -t TABLE, --table TABLE
                        color table (optional)
  -d DEFAULT_COLORS, --default_colors DEFAULT_COLORS
                        TSV file that specifies default color Configurator (optional; default: data/default_colors.tsv)
  -n NT, --nt NT        dinucleotide (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
  --species SPECIES     Species name (optional; e.g. "<i>Escherichia coli</i>", "<i>Ca.</i> Hepatoplasma crinochetorum")
  --strain STRAIN       Strain/isolate name (optional; e.g. "K-12", "Av")
  -k FEATURES, --features FEATURES
                        Comma-separated list of feature keys to draw (default: CDS,tRNA,rRNA,repeat_region)
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "black")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (float; default: 0)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (str; default: "gray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (float; default: 1.0)
  -f FORMAT, --format FORMAT
                        Comma-separated list of output file formats (default: png)
  --suppress_gc         Suppress GC content track (default: False).
  --suppress_skew       Suppress GC skew track (default: False).
```
#### <i>Haemophilus influenzae</i>
`gbdraw` automatically identifies and displays the organism and strain name from the sequence record. However, these names are not italicized by default. For example:
```bash
gbdraw circular -i GCF_000931575.1_ASM93157v1_genomic.gbff -o Haemophilus_influenzae
```
![hinfluenzae](https://github.com/satoshikawato/gbdraw/blob/main/examples/Haemophilus_influenzae.png)
#### <i>Escherichia coli</i> K-12
To italicize a portion of the organism name, you can use the <i></i> tags in the --species and --strain parameters. This will format the specified text in italics. The following command will render the species name "_Escherichia coli_" in italics, while keeping "K-12" in standard text (the organim name will be overridden):
```bash
gbdraw circular -i NC_000913.gb --species "<i>Escherichia coli</i>" --strain "K-12"
```
![ecoli](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000913.png)

#### <i>Vibrio cholerae</i> Strain O395 (GCF_000016245.1)
For GenBank files containing multiple entries, `gbdraw` saves each entry as a separate file. Here's how to do it for the _Vibrio cholerae_ strain O395 genome, which has two chromosomes:
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/016/245/GCF_000016245.1_ASM1624v1/GCF_000016245.1_ASM1624v1_genomic.gbff.gz # download genome
gunzip GCF_000016245.1_ASM1624v1_genomic.gbff.gz # extract file
gbdraw circular -i GCF_000016245.1_ASM1624v1_genomic.gbff --species "<i>Vibrio cholerae</i>" --strain "O395" -f svg  # Draw genome; results in "NC_009457.svg" for Chromosome I and "NC_009456.svg" for Chromosome II
```
<i>Vibrio cholerae</i> Chromosome I
![Vibrio cholerae chromosome I](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009457.png)
<i>Vibrio cholerae</i> Chromosome II
![Vibrio cholerae chromosome II](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_009456.png)
### Linear genome
`gbdraw linear`
```bash
$ gbdraw linear -h
usage: gbdraw [-h] -i [INPUT ...] [-b [BLAST ...]] [-t TABLE] [-d DEFAULT_COLORS] [-o OUTPUT] [-n NT] [-w WINDOW] [-s STEP]
              [--separate_strands] [--show_gc] [--align_center] [--evalue EVALUE] [--bitscore BITSCORE] [--identity IDENTITY]
              [-k FEATURES] [--block_stroke_color BLOCK_STROKE_COLOR] [--block_stroke_width BLOCK_STROKE_WIDTH]
              [--line_stroke_color LINE_STROKE_COLOR] [--line_stroke_width LINE_STROKE_WIDTH] [-f FORMAT]

Generate plot in PNG/PDF/SVG/PS/EPS.

options:
  -h, --help            show this help message and exit
  -i [INPUT ...], --input [INPUT ...]
                        genbank (required)
  -b [BLAST ...], --blast [BLAST ...]
                        input BLAST result file in tab-separated format (-outfmt 6 or 7) (optional)
  -t TABLE, --table TABLE
                        color table (optional)
  -d DEFAULT_COLORS, --default_colors DEFAULT_COLORS
                        TSV file that specifies default color Configurator (optional; default: data/default_colors.tsv)
  -o OUTPUT, --output OUTPUT
                        output file prefix (default: out)
  -n NT, --nt NT        dinucleotide skew (default: GC).
  -w WINDOW, --window WINDOW
                        window size (default: 1000)
  -s STEP, --step STEP  step size (default: 100)
  --separate_strands    separate forward and reverse strands (default: False). Features of undefined strands are shown on the
                        forward strand.
  --show_gc             plot GC content below genome (default: False).
  --align_center        Align genomes to the center (default: False).
  --evalue EVALUE       evalue threshold (default=1e-2)
  --bitscore BITSCORE   bitscore threshold (default=50)
  --identity IDENTITY   identity threshold (default=0)
  -k FEATURES, --features FEATURES
                        Comma-separated list of feature keys to draw (default: CDS,tRNA,rRNA,repeat_region)
  --block_stroke_color BLOCK_STROKE_COLOR
                        Block stroke color (str; default: "black")
  --block_stroke_width BLOCK_STROKE_WIDTH
                        Block stroke width (float; default: 0)
  --line_stroke_color LINE_STROKE_COLOR
                        Line stroke color (str; default: "gray")
  --line_stroke_width LINE_STROKE_WIDTH
                        Line stroke width (float; default: 1.0)
  -f FORMAT, --format FORMAT
                        Comma-separated list of output file formats (default: png)
```
### Human herpesvirus 6 (HHV-6)
`gbdraw` can draw linear genomes and pairwise matches depicting similar genomic regions.
Perform a pariwise BLASTN analysis between the genomes of HHV-6A and HHV-6B:
```bash
blastn -query NC_000898.fasta -subject NC_001664.fasta -outfmt 7 -out NC_000898_NC_001664.blastn.out
```
Use `gbdraw linear` to visualize the comparison, aligning genomes to the center and separating forward and reverse strands:
```bash
gbdraw linear -i NC_000898.gb NC_001664.gb -b NC_000898_NC_001664.blastn.out --align_center --separate_strands -o HHV-6
```
![HHV-6](https://github.com/satoshikawato/gbdraw/blob/main/examples/HHV-6.png)

### <i>Saccharomyces cerevisiae</i> (budding yeast)
Although primarily focused on smaller genomes, `gbdraw` <i>can</i> handle eukaryotic genomes like the budding yeast, <i>S. cerevisiae</i>:
```bash
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64/GCF_000146045.2_R64_genomic.gbff.gz
gunzip GCF_000146045.2_R64_genomic.gbff.gz
gbdraw linear -i GCF_000146045.2_R64_genomic.gbff --separate_strands --show_gc -o yeast # outputs yeast.png
```
![yeast](https://github.com/satoshikawato/gbdraw/blob/main/examples/yeast.png)
## Advanced customization
### Customizing colors with configuration files
`gbdraw` allows users to personalize feature colors. This is done by providing a custom tab-separated file that specifies the desired color codes. Additionally, the colors and widths of the strokes can be fine-tuned using command-line arguments like `--block_stroke_width`.
Colors can be specified using either of the following methods:
Color Names: You can use any of the [147 color names defined by the SVG specification](https://johndecember.com/html/spec/colorsvg.html).
Hexadecimal Format: For more precise color control, you can express colors in [hexadecimal format](https://htmlcolorcodes.com/).
#### Default color table
The default color table (`/site-packages/gbdraw/data/default_colors.txt`) is as follows:
```default_colors.txt
CDS     #47b8f8
rRNA    #009e73
tRNA    #e69f00
repeat_region   #d3d3d3
misc_feature    #d3d3d3
default #d3d3d3
```
The following `modified_default_colors.txt`
```modified_default_colors.txt
CDS	#993366
rRNA	#00ff00
tRNA	#228b22
repeat_region	#a6acb3
```
#### Feature-specific color table
For more precise control, especially for CDS features, you can specify colors based on specific attributes like `product` or `note`.
```feature_specific_color_table.txt
CDS	product	hypothetical protein	#b3b3b3
CDS	note	possible pseudo*	#dba2bf
```
Here, any CDS feature with unknown function (`hypothetical protein`) is colored gray. Possible pseudogenes are colored light purple.
```bash
gbdraw circular -i AP027131.gb --species "<i>Ca.</i> Hepatoplasma vulgare" --strain "Av-JP" -d modified_default_colors.txt -t feature_specific_color_table.txt -o Av-JP -f svg # PNG and PDF outputs are bugged;"Ca." cannot be placed correctly
```
![Av-JP](https://github.com/satoshikawato/gbdraw/blob/main/examples/Av-JP.svg)
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
MjeNMV.gb \
MelaMJNV.gb \
PemoMJNVA.gb \
PeseMJNV.gb \
PemoMJNVB.gb \
LvMJNV.gb \
TrcuMJNV.gb \
MellatMJNV.gb \
MeenMJNV.gb \
MejoMJNV.gb \
-b \
MjeNMV.MelaMJNV.tblastx.out \
MelaMJNV.PemoMJNVA.tblastx.out \
PemoMJNVA.PeseMJNV.tblastx.out \
PeseMJNV.PemoMJNVB.tblastx.out \
PemoMJNVB.LvMJNV.tblastx.out \
LvMJNV.TrcuMJNV.tblastx.out \
TrcuMJNV.MellatMJNV.tblastx.out \
MellatMJNV.MeenMJNV.tblastx.out \
MeenMJNV.MejoMJNV.tblastx.out \
-t color_table.txt \
-d modified_default_colors.tsv \
--align_center \
--separate_strands \
-o majani 
```
![majaniviruses](https://github.com/satoshikawato/gbdraw/blob/main/examples/majani.svg)

## Planned features
- Feature color code legends (will be added in the next update)
- Dynamic scaling of the text (will be added in the next update)
- Feature labels (hopefully)
- Multiple tracks to visualize overlapping features (overlapping genes, transcript isoforms etc.)

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

The core functionality of gbdraw has evolved from [a set of Python scripts](https://github.com/satoshikawato/bio_small_scripts/) written back in 2022. Key among these are:

[plot_circular_genome.py](https://github.com/satoshikawato/bio_small_scripts/blob/main/plot_circular_genome.py)

[plot_linear_genome.py](https://github.com/satoshikawato/bio_small_scripts/blob/main/plot_linear_genome.py)

