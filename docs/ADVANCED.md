[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)
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
gbdraw circular --gbk LC738868.gb -o LC738868_middle_separate_strands -f svg,png --block_stroke_width 1 --block_stroke_color gray --show_labels --track_type middle --separate_strands -t custom_color_table.tsv -d modified_default_colors.tsv 
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
--gbk \
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
gbdraw circular --gbk AP027280.gb -f svg --block_stroke_width 1 --block_stroke_color gray --track_type spreadout --show_labels
```
![WSSV](https://github.com/satoshikawato/gbdraw/blob/main/examples/AP027280.svg)
```bash
gbdraw circular --gbk NC_012920.gb -f svg --block_stroke_width 2 --block_stroke_color gray  --show_labels -w 100 -s 10
```
![HsmtDMA](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_012920.svg)

#### Label size adjustment
The font size of the Label text can be adjusted via `--label_font_size` option.
#### Label text blacklist 
You can selectively hide some labels via `--label_blacklist` option. For example, supressing `hypothetical protein` and other less-informative labels reduces label redundancy and makes it easier to visualize larger genomes without layout issues.
`--label_blacklist` accepts comma-separated keywords or path to a file containing single keyword per line.
```bash
gbdraw circular --gbk AP027280.gb \
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
--gbk AP027078.gb \
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

#### Lable whitelist

```bash
gbdraw circular \
--gbk NC_010162.gb \
-f svg \
--palette edelweiss \
--show_labels \
--separate_strands \
-t NC_010162.feature-specific_table.tsv \
--label_whitelist NC_010162.whitelist.tsv \
-t NC_010162.feature-specific_table.tsv \
--species "<i>Sorangium cellulosum</i>" \
--strain "So ce56"
```
##### Lable whitelist example
| feature type | target qualifier | qualifier value regex (Python) |
| ------ | ------- | ------- |
| CDS | old_locus_tag | sce4138 |
| CDS | old_locus_tag | sce4137 |

[NC_010162.whitelist.tsv](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162.whitelist.tsv) (part)
```
CDS	old_locus_tag	sce4138
CDS	old_locus_tag	sce4137
CDS	old_locus_tag	sce4136
CDS	old_locus_tag	sce4135
CDS	old_locus_tag	sce4134
CDS	old_locus_tag	sce4133
CDS	old_locus_tag	sce4132
...
```
[NC_010162.feature-specific_table.tsv](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162.feature-specific_table.tsv) (part)
```
CDS	old_locus_tag	sce4138	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4137	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4136	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4135	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4134	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4133	#a4d8a7	Chivosazol biosynthesis
CDS	old_locus_tag	sce4132	#a4d8a7	Chivosazol biosynthesis
...
```
![NC_010162_edelweiss](https://github.com/satoshikawato/gbdraw/blob/main/examples/NC_010162_edelweiss.svg)


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
--gbk NC_012920.gb \
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

[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)