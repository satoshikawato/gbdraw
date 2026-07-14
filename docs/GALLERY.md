[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | **Gallery** | [FAQ](./FAQ.md) | [About](./ABOUT.md)
# Gallery
This gallery contains representative `gbdraw` outputs and the commands used to generate them. For shorter copy-paste examples, see [Recipes](./RECIPES.md).

Interactive SVG versions of selected examples are available at [https://gbdraw.app/gallery/](https://gbdraw.app/gallery/). GitHub previews static SVG images, so use the public gallery when you want to inspect JavaScript-enabled SVG output. Selected gallery entries also include a Tutorial tab for reproducing the same figure in the web app.

The interactive gallery assets are hosted online and are not bundled with local `gbdraw` installs. Running `gbdraw gui` locally still keeps genome analysis on your own machine.

For circular examples, `--track_type` names the circular layout preset. Custom track slots use explicit geometry instead, and the circular axis remains fixed.

When converting long gallery commands into TSV manifests, see the `--records_table`, `--conservation_table`, and `--circular_track_table` examples in the [CLI Reference](./CLI_Reference.md#tsv-manifest-inputs).

#### White spot syndrome virus genome comparison

Interactive SVG session output: [https://gbdraw.app/gallery/#WSSV_genome_comparison](https://gbdraw.app/gallery/#WSSV_genome_comparison)

Web tutorial: open the [WSSV gallery entry](https://gbdraw.app/gallery/#WSSV_genome_comparison) and select the Tutorial tab.

The tutorial builds a circular WSSV conservation-ring plot in the browser: it runs LOSAT `blastn`, uploads comparison FASTA files, sets ring labels and colors, checks feature and legend popups, and exports a session.


#### <i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018 complete genome

Interactive SVG session output: [https://gbdraw.app/gallery/#Vnig_TUMSAT-TG-2018](https://gbdraw.app/gallery/#Vnig_TUMSAT-TG-2018)

Web tutorial: open the [Vibrio multi-record gallery entry](https://gbdraw.app/gallery/#Vnig_TUMSAT-TG-2018) and select the Tutorial tab.

This circular multi-record example starts from a hosted gallery session and interactive SVG output.

The tutorial uses one GBFF input with six records, sets up the multi-record canvas, assigns rows, checks feature popups, and exports a session.

#### Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.

Interactive SVG session output: [https://gbdraw.app/gallery/#BGC0000708-BGC0000713](https://gbdraw.app/gallery/#BGC0000708-BGC0000713)

Web tutorial: open the [BGC gallery entry](https://gbdraw.app/gallery/#BGC0000708-BGC0000713) and select the Tutorial tab.

The tutorial runs a browser LOSATP orthogroup search, colors features with antiSMASH `gene_kind` rules, rotates first-record labels, checks feature and orthogroup popups, edits the legend, and exports a session.

#### Majaniviruses multi-genome comparison

Interactive SVG session output: [https://gbdraw.app/gallery/#majanivirus_orthogroup](https://gbdraw.app/gallery/#majanivirus_orthogroup)

Web tutorial: open the [majanivirus gallery entry](https://gbdraw.app/gallery/#majanivirus_orthogroup) and select the Tutorial tab.

The tutorial runs a nine-record LOSATP orthogroup comparison in the browser, applies product-based color rules, sets record labels, checks match popups, and exports a session.

<details><summary>Expand to see the script</summary>

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
--gbk \
./in_gbk/MjeNMV.gb \
./in_gbk/MelaMJNV.gb \
./in_gbk/PemoMJNVA.gb \
./in_gbk/PeseMJNV.gb \
./in_gbk/PemoMJNVB.gb \
./in_gbk/LvMJNV.gb \
./in_gbk/TrcuMJNV.gb \
./in_gbk/MellatMJNV.gb \
./in_gbk/MeenMJNV.gb \
./in_gbk/MejoMJNV.gb \
-b \
./in_fna/MjeNMV.MelaMJNV.tblastx.out \
./in_fna/MelaMJNV.PemoMJNVA.tblastx.out \
./in_fna/PemoMJNVA.PeseMJNV.tblastx.out \
./in_fna/PeseMJNV.PemoMJNVB.tblastx.out \
./in_fna/PemoMJNVB.LvMJNV.tblastx.out \
./in_fna/LvMJNV.TrcuMJNV.tblastx.out \
./in_fna/TrcuMJNV.MellatMJNV.tblastx.out \
./in_fna/MellatMJNV.MeenMJNV.tblastx.out \
./in_fna/MeenMJNV.MejoMJNV.tblastx.out \
-t majani_custom_color_table.tsv \
-d modified_default_colors.tsv \
--block_stroke_width 1 \
--block_stroke_color gray \
--align_center \
--separate_strands \
-o majani -f svg
```

</details>

![majaniviruses](https://github.com/satoshikawato/gbdraw/blob/main/examples/majani.svg)

#### Hepatoplasmataceae five-genome comparison

Interactive SVG session outputs: [collinear](https://gbdraw.app/gallery/#hepatoplasmataceae_collinear), [orthogroup](https://gbdraw.app/gallery/#hepatoplasmataceae_orthogroup)

Web tutorials: open the [collinear gallery entry](https://gbdraw.app/gallery/#hepatoplasmataceae_collinear) or [orthogroup gallery entry](https://gbdraw.app/gallery/#hepatoplasmataceae_orthogroup) and select the Tutorial tab.

<details><summary>Expand to see the script</summary>

```bash
gbdraw linear \
--gbk \
AP027078.gb \
AP027131.gb \
AP027133.gb \
AP027132.gb \
NZ_CP006932.gb \
-b \
AP027078_AP027131.tblastx.out \
AP027131_AP027133.tblastx.out \
AP027133_AP027132.tblastx.out \
AP027132_NZ_CP006932.tblastx.out \
--align_center \
--separate_strands \
--show_gc \
--show_skew \
--block_stroke_width 1 \
--block_stroke_color gray \
--palette default \
-f svg \
-o hepatoplasmataceae_default
```

</details>

![hepatoplasmataceae_default.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/hepatoplasmataceae_default.svg)

#### Human mitochondrial genome (feature qualifier priority)

Interactive SVG session output: [https://gbdraw.app/gallery/#HmmtDNA_ATskew](https://gbdraw.app/gallery/#HmmtDNA_ATskew)

Web tutorial: open the [human mitochondrial AT skew gallery entry](https://gbdraw.app/gallery/#HmmtDNA_ATskew) and select the Tutorial tab.

The tutorial uses custom track slots to add an AT-skew ring, switches CDS labels to gene names, checks feature popups, edits the legend, and exports a session.

<details><summary>Expand to see the script</summary>

The `qualifier_priority.tsv` file passed with `--qualifier_priority` chooses the label qualifier for each feature type. Here, CDS labels use `gene`; other feature types keep the default qualifier order:

```qualifier_priority.tsv
CDS gene
```
| feature type | qualifier |
| ------ | ------- |
| CDS | gene |

```bash
gbdraw circular \
--gbk NC_012920.gb \
-f svg --track_type middle \
--species "<i>Homo sapiens</i>" \
--block_stroke_width 2 \
--axis_stroke_width 5 \
--labels both \
--qualifier_priority qualifier_priority.tsv \
-o NC_012920_middle_qualifier_priority_inner_axis5_def28_italic \
--definition_font_size 28
```
</details>

![NC_012920_middle_qualifier_priority_inner_axis5_def28_italic](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_012920_middle_qualifier_priority_inner_axis5_def28_italic.svg)

#### <i>Nicotiana tabacum</i> chloroplast genome

<details><summary>Expand to see the script</summary>


```bash
gbdraw circular \
--gbk NC_001879.gbk \
--separate_strands \
-f svg \
-o NC_001879_color \
-k CDS,rRNA,tRNA,tmRNA,ncRNA,misc_RNA,rep_origin \
-t chloroplast_specific_table.tsv \
--block_stroke_width 1 \
--block_stroke_color black \
--axis_stroke_width 3 \
--line_stroke_width 2 \
--suppress_gc \
--suppress_skew \
-p default \
--track_type tuckin \
--labels both \
--qualifier_priority qualifier_priority.tsv \
--outer_label_x_radius_offset 0.90 \
--outer_label_y_radius_offset 0.90 \
--inner_label_x_radius_offset 0.975 \
--inner_label_y_radius_offset 0.975 \
--species "<i>Nicotiana tabacum</i>" \
--definition_font_size 28 \
--legend upper_left
```
</details>

![NC_001879_color.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_001879_color.svg)

#### Lambda phage

<details><summary>Expand to see the script</summary>


```bash
gbdraw linear \
--gbk NC_001416.gb \
-o NC_001416 \
-f svg \
--show_labels all \
--separate_strands \
--legend left \
-d cds_white.tsv \
-t lambda_specific_table.tsv \
--block_stroke_width 2 \
--axis_stroke_width 5 \
--definition_font_size 24
```

</details>

![NC_001416.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_001416.svg)


#### <i>Ca.</i> Sukunaarchaeum mirabile M16-5

<details><summary>Expand to see the script</summary>


```bash
gbdraw circular \
--gbk M16-5.gb \
--separate_strands \
-f svg \
-o M16-5_fugaku \
--block_stroke_width 1 \
--axis_stroke_width 1 \
-p fugaku \
--track_type middle \
--species "<i>Ca.</i> Sukunaarchaeum mirabile" \
--definition_font_size 22 \
--legend upper_right
```

</details>

![M16-5_fugaku.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/M16-5_fugaku.svg)


#### <i>Pandoravirus salinus</i>

<details><summary>Expand to see the script</summary>

```bash
gbdraw circular \
--gbk Pandoravirus_salinus.gb \
--separate_strands \
-f svg \
-o Pandoravirus_salinus_forest \
-p forest \
--track_type tuckin \
--species "<i>Pandoravirus salinus</i>" \
--definition_font_size 22 \
--legend upper_right
```

</details>

![Pandoravirus_salinus_forest.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/Pandoravirus_salinus_forest.svg)

#### <i>Candidatus</i> Pelagibacter ubique HTCC1062

<details><summary>Expand to see the script</summary>

```bash
gbdraw circular \
--gbk NC_007205.gb \
-f svg \
--separate_strands \
--species "<i>Ca. </i> Pelagibacter ubique" \
--strain "HTCC1062" \
--legend none \
--palette oceanic_voyage \
-o NC_007205_oceanic_voyage
```

</details>

![NC_007205_oceanic_voyage.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_007205_oceanic_voyage.svg)

#### <i>Prochlorococcus marinus</i> CCMP1375

<details><summary>Expand to see the script</summary>


```bash
gbdraw circular \
--gbk NC_005042.gb \
-f svg \
--separate_strands \
--species "<i>Prochlorococcus marinus</i>" \
--strain "CCMP1375" \
--legend none \
--palette pine_reflection \
-o NC_005042_pine_reflection
```

</details>

![NC_005042_pine_reflection.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_005042_pine_reflection.svg)

#### <i>Flavobacterium columnare</i> ATCC 49512

<details><summary>Expand to see the script</summary>


```bash
gbdraw circular \
--gbk NC_016510.gb \
-f svg \
--separate_strands \
--species "<i>Flavobacterium columnare</i>" \
--strain "ATCC 49512" \
--legend none \
--palette mint \
-o NC_016510_mint
```

</details>

![NC_016510_mint.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_016510_mint.svg)

#### <i>Thermus aquaticus</i> Y51MC23

<details><summary>Expand to see the script</summary>

```bash
gbdraw circular \
--gbk NZ_CP010822.gb \
-f svg \
--separate_strands \
--species "<i>Thermus aquaticus</i>" \
--strain "Y51MC23" \
--legend none \
--palette orange \
-o NZ_CP010822_orange
```
</details>

![NZ_CP010822_orange.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NZ_CP010822_orange.svg)

#### <i>Helicobacter pylori</i> J99

<details><summary>Expand to see the script</summary>

```bash
gbdraw circular \
--gbk NC_000921.gb \
-f svg \
--separate_strands \
--species "<i>Helicobacter pylori</i>" \
--strain "J99" \
--legend none \
--palette spring \
-o NC_000921_spring
```

</details>

![NC_000921_spring.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000921_spring.svg)


#### <i>Mycobacterium tuberculosis</i> H37Rv

<details><summary>Expand to see the script</summary>

```bash
gbdraw circular \
--gbk NC_000962.gb \
-f svg \
--separate_strands \
--species "<i>Mycobacterium tuberculosis</i>" \
--strain "H37Rv" \
--legend none \
--palette psyche \
-o NC_000962_psyche
```

</details>

![NC_000962_psyche.svg](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_000962_psyche.svg)


[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Recipes](./RECIPES.md) | [CLI Reference](./CLI_Reference.md) | **Gallery** | [FAQ](./FAQ.md) | [About](./ABOUT.md)
