# Internal CLI and Python scenario evidence

This registry stores literal executable blocks and runner-checked TSV inputs for
internal H-CLI and H-PY regression scenarios. It is not public user
documentation. Public procedures and contracts live in Tutorials and Technical
documentation.

## Runner-owned TSV literals
```tsv
gbk	record_label	record_subtitle	record_id	order	row	column
../HmmtDNA.gbk	Human mitochondrion	Complete RefSeq record	NC_012920.1	1	1	1
../NC_002333.2.gb	Zebrafish mitochondrion	Complete RefSeq record	NC_002333.2	2	1	2
../NC_024511.2.gb	Fruit fly mitochondrion	Complete RefSeq record	NC_024511.2	3	2	1
../NC_001328.1.gb	Nematode mitochondrion	Complete RefSeq record	NC_001328.1	4	2	2
```

```tsv
blast	query	subject
../lambda-de3.losatn.tsv	NC_001416.1	NC_042057.1
```

```tsv
blast	label	color
../danio-human.tlosatx.tsv	Danio rerio (NC_002333.2)	#4E79A7
../drosophila-human.tlosatx.tsv	Drosophila melanogaster (NC_024511.2)	#F28E2B
../caenorhabditis-human.tlosatx.tsv	Caenorhabditis elegans (NC_001328.1)	#59A14F
```

```tsv
id	renderer	side	r	w	inner_gap_px	outer_gap_px	z	params
ticks	ticks	outside						tick_label_layout=label_out_tick_in
features	features	axis						
gc_content	dinucleotide_content	inside		0.10	3	3		nt=GC,legend_label=GC content
gc_skew	dinucleotide_skew	inside		0.10	3	3		nt=GC,legend_label=GC skew
at_skew	dinucleotide_skew	inside		0.10	3	3		nt=AT,positive_color=#deaf6e,negative_color=#7294e3,legend_label=AT skew
```

```tsv
CDS	gene	^ND(4L|[1-6])$	#3B82F6	NADH dehydrogenase
CDS	gene	^COX[1-3]$	#EF4444	Cytochrome c oxidase
CDS	gene	^ATP[68]$	#F59E0B	ATP synthase
CDS	gene	^CYTB$	#8B5CF6	Cytochrome b
rRNA	gene	^RNR[12]$	#10B981	Ribosomal RNA
```

```tsv
CDS	gene	^(ND1|COX2|ATP8|ATP6|COX3|CYTB)$
rRNA	gene	^RNR[12]$
```

```tsv
record_id	feature_type	qualifier	value	label_text
NC_012920.1	CDS	label	^ND1$	Complex I (ND1)
NC_012920.1	CDS	label	^COX2$	Oxidase II
NC_012920.1	rRNA	label	^s-rRNA$	12S rRNA
NC_012920.1	rRNA	label	^l-rRNA$	16S rRNA
```

## Executable blocks

<!-- executable:H-CLI-09:start -->
```bash
gbdraw circular \
  --gbk AP027133.gb \
  --depth_track AP027133.DRR394922.depth-1kb.tsv \
  --depth_track_label 'DRR394922 mean depth' \
  --depth_track_color '#2563EB' \
  --depth_window 1 \
  --depth_step 1000 \
  --depth_min 0 \
  --depth_max 80 \
  --no_depth_log_scale \
  --show_depth_axis \
  --show_depth_ticks \
  --depth_large_tick_interval 20 \
  --depth_small_tick_interval 10 \
  --gc \
  --skew \
  --window 1000 \
  --step 1000 \
  --gc_content_mode percent \
  --gc_content_min_percent 10 \
  --gc_content_max_percent 55 \
  --show_gc_content_axis \
  --show_gc_content_ticks \
  --gc_content_large_tick_interval 10 \
  --gc_content_small_tick_interval 5 \
  --circular_track_slot 'ticks:ticks@side=outside,tick_label_layout=label_out_tick_in' \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split,w=48px' \
  --circular_track_slot 'depth_1:depth@track_index=0,side=inside,w=52px,legend_label=DRR394922 mean depth' \
  --circular_track_slot 'gc_content:dinucleotide_content@side=inside,w=42px,nt=GC,legend_label=GC content (%)' \
  --circular_track_slot 'gc_skew:dinucleotide_skew@side=inside,w=34px,nt=GC,legend_label=GC skew' \
  --circular_track_slot 'at_skew:dinucleotide_skew@side=inside,w=34px,nt=AT,positive_color=#DEAF6E,negative_color=#7294E3,legend_label=AT skew' \
  --separate_strands \
  --labels none \
  --plot_title 'Depth, GC, and skew tracks' \
  --plot_title_position top \
  --plot_title_font_size 20 \
  --legend right \
  -o cli_quantitative_tracks \
  -f svg
```
<!-- executable:H-CLI-09:end -->

<!-- executable:H-CLI-10:start -->
```bash
gbdraw circular \
  --gbk NC_001879.gbk \
  --annotation_table nicotiana-tabacum-regions.tsv \
  -d modified_default_colors.tsv \
  --qualifier_priority qualifier_priority.tsv \
  --circular_track_slot 'ticks:ticks@side=outside,tick_label_layout=label_out_tick_in' \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split,w=48px' \
  --circular_track_slot 'plastome_regions:annotations@set_id=plastome_regions,side=inside,w=30px,show_labels=true,padding_px=1,overflow=compress' \
  --circular_track_slot 'gc_content:dinucleotide_content@side=inside,w=36px,nt=GC,legend_label=GC content' \
  --circular_track_slot 'gc_skew:dinucleotide_skew@side=inside,w=32px,nt=GC,legend_label=GC skew' \
  --separate_strands \
  --track_type tuckin \
  --labels none \
  --plot_title 'Nicotiana tabacum plastome structure' \
  --plot_title_position bottom \
  --legend right \
  -o cli_annotations_slots \
  -f svg
```
<!-- executable:H-CLI-10:end -->

<!-- executable:H-CLI-04:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb BGC0000708.gbk BGC0000713.gbk \
  --record_id NC_001416.1 \
  --record_id BGC0000708 \
  --record_id BGC0000713 \
  --region NC_001416.1:5001-35500 \
  --reverse_complement false \
  --reverse_complement false \
  --reverse_complement true \
  --multi_record_position '#1@1' \
  --multi_record_position '#2@2' \
  --multi_record_position '#3@2' \
  --record_label 'Lambda selected region' \
  --record_label 'Lividomycin cluster' \
  --record_label 'Ribostamycin cluster' \
  --record_subtitle 'NC_001416.1 positions 5,001–35,500' \
  --record_subtitle 'Complete BGC0000708 record' \
  --record_subtitle 'Complete BGC0000713 reverse complement' \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --show_labels all \
  --track_layout above \
  --scale_style ruler \
  --ruler_on_axis \
  --linear_record_gap 28 \
  --keep_definition_left_aligned \
  --legend right \
  -o linear_layout_cli \
  -f svg
```
<!-- executable:H-CLI-04:end -->

<!-- executable:H-CLI-03:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk NC_002333.2.gb NC_024511.2.gb NC_001328.1.gb \
  --multi_record_canvas \
  --multi_record_size_mode equal \
  --multi_record_position '#1@1' \
  --multi_record_position '#2@1' \
  --multi_record_position '#3@2' \
  --multi_record_position '#4@2' \
  --circular_track_order ticks,features,gc_content,gc_skew \
  --circular_track_axis_index 1 \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels out \
  --label_font_size 10 \
  --plot_title 'Four complete metazoan mitochondrial genomes' \
  --plot_title_position bottom \
  --legend right \
  -o multi_record_circular_cli \
  -f svg
```
<!-- executable:H-CLI-03:end -->

<!-- executable:H-CLI-07:start -->
```bash
gbdraw linear \
  --gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk \
  --protein_blastp_mode orthogroup \
  --losatp_threads 1 \
  --identity 30 \
  --align_orthogroup_feature CAG38695.1 \
  -k CDS,rRNA,tRNA,tmRNA,ncRNA,repeat_region \
  -p orange \
  -d BGC0000708-BGC0000713_default_colors.tsv \
  -t BGC0000708-BGC0000713_specific_colors.tsv \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --show_labels first \
  --label_placement above_feature \
  --label_rotation 45 \
  --pairwise_match_style curve \
  --scale_style ruler \
  --plot_title 'Aminoglycoside biosynthetic gene clusters from <i>Streptomyces</i> spp.' \
  --keep_definition_left_aligned \
  --block_stroke_width 2 \
  --block_stroke_color '#262626' \
  --line_stroke_width 2 \
  --axis_stroke_width 5 \
  --legend_box_size 20 \
  --legend_font_size 20 \
  --label_font_size 18 \
  --feature_height 75 \
  --ruler_label_font_size 20 \
  --definition_line_style name:size=20,weight=bold \
  --definition_line_style subtitle:size=20 \
  --definition_line_style 'accession:size=20,color=#7b7c7d' \
  --definition_line_style 'length:size=20,color=#7b7c7d' \
  --legend bottom \
  -o cli_losatp_groups \
  -f svg
```
<!-- executable:H-CLI-07:end -->

<!-- executable:H-CLI-08:start -->
```bash
gbdraw linear \
  --gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk \
  --protein_blastp_mode collinear \
  --losatp_threads 1 \
  --identity 30 \
  --collinear_search_scope adjacent \
  --collinear_min_anchors 2 \
  --collinear_max_unit_gap 1 \
  --collinear_max_diagonal_drift 1 \
  --collinear_color_mode orientation_identity \
  -k CDS,rRNA,tRNA,tmRNA,ncRNA,repeat_region \
  -p orange \
  -d BGC0000708-BGC0000713_default_colors.tsv \
  -t BGC0000708-BGC0000713_specific_colors.tsv \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --show_labels first \
  --label_placement above_feature \
  --label_rotation 45 \
  --pairwise_match_style curve \
  --scale_style ruler \
  --plot_title 'Collinear protein-match blocks across five BGC records' \
  --keep_definition_left_aligned \
  --block_stroke_width 2 \
  --block_stroke_color '#262626' \
  --line_stroke_width 2 \
  --axis_stroke_width 5 \
  --legend_box_size 20 \
  --legend_font_size 20 \
  --label_font_size 18 \
  --feature_height 75 \
  --ruler_label_font_size 20 \
  --definition_line_style name:size=20,weight=bold \
  --definition_line_style subtitle:size=20 \
  --definition_line_style 'accession:size=20,color=#7b7c7d' \
  --definition_line_style 'length:size=20,color=#7b7c7d' \
  --legend bottom \
  -o cli_losatp_collinear \
  -f svg
```
<!-- executable:H-CLI-08:end -->

<!-- executable:H-CLI-05:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb NC_042057.1.gb \
  --record_id NC_001416.1 \
  --record_id NC_042057.1 \
  --blast lambda-de3.losatn.tsv \
  --bitscore 50 \
  --evalue 1e-2 \
  --identity 0 \
  --alignment_length 0 \
  --pairwise_match_style curve \
  --show_labels none \
  --track_layout above \
  --scale_style ruler \
  --ruler_on_axis \
  --legend right \
  -o linear_precomputed_comparison \
  -f svg
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --conservation_blast danio-human.tlosatx.tsv drosophila-human.tlosatx.tsv caenorhabditis-human.tlosatx.tsv \
  --conservation_reference subject \
  --conservation_labels 'Danio rerio (NC_002333.2)' 'Drosophila melanogaster (NC_024511.2)' 'Caenorhabditis elegans (NC_001328.1)' \
  --conservation_colors '#4E79A7' '#F28E2B' '#59A14F' \
  --bitscore 50 \
  --evalue 1e-5 \
  --identity 40 \
  --alignment_length 50 \
  --conservation_ring_width 18 \
  --conservation_ring_gap 4 \
  --track_type middle \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels out \
  --plot_title 'Complete metazoan mitochondrial TLOSATX evidence' \
  --plot_title_position top \
  --legend right \
  -o circular_conservation_ring \
  -f svg
```
<!-- executable:H-CLI-05:end -->

<!-- executable:H-CLI-13:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --species '<i>Homo sapiens</i>' \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels both \
  --track_type middle \
  --plot_title 'Human mitochondrial genome export set' \
  --plot_title_position top \
  --plot_title_font_size 20 \
  --legend right \
  -o cli_export \
  -f svg,interactive_svg,png,pdf,eps,ps
```
<!-- executable:H-CLI-13:end -->

<!-- executable:H-CLI-06:start -->
```bash
gbdraw linear \
  --gbk BGC0000708.gbk BGC0000709.gbk BGC0000711.gbk BGC0000712.gbk BGC0000713.gbk \
  --protein_blastp_mode pairwise \
  --losatp_threads 1 \
  --protein_blastp_max_hits 1 \
  --protein_blastp_output cli_losatp_pairwise.tsv \
  --identity 30 \
  --pairwise_match_style curve \
  --show_labels none \
  --track_layout above \
  --scale_style ruler \
  --ruler_on_axis \
  --align_center \
  --plot_title 'LOSATP Pairwise protein matches across five BGC records' \
  --legend bottom \
  -o cli_losatp_pairwise \
  -f svg
```
<!-- executable:H-CLI-06:end -->

<!-- executable:H-CLI-12:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  --species '<i>Homo sapiens</i>' \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels both \
  --track_type middle \
  --session_output cli_session.json \
  -o cli_session_roundtrip \
  -f svg

gbdraw circular \
  --session cli_session.json \
  --session_output cli_session.json.gz \
  --overwrite \
  -o cli_session_roundtrip \
  -f svg
```
<!-- executable:H-CLI-12:end -->

<!-- executable:H-CLI-11:start -->
```bash
gbdraw circular \
  --gbk HmmtDNA.gbk \
  -k CDS,rRNA,tRNA \
  --palette colorblind \
  --table tables/presentation_colors.tsv \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --label_whitelist tables/presentation_labels.tsv \
  --label_table tables/presentation_label_overrides.tsv \
  --feature_visibility_table HmmtDNA_feature_visibility.tsv \
  --feature_shape CDS=arrow \
  --feature_shape rRNA=rectangle \
  --feature_shape tRNA=arrow \
  --feature_shape D-loop=underlay \
  --arrow_head_length_ratio auto \
  --arrow_shaft_width_ratio 0.72 \
  --resolve_overlaps \
  --separate_strands \
  --track_type spreadout \
  --labels both \
  --label_rendering auto \
  --block_stroke_color '#1F2937' \
  --block_stroke_width 1.5 \
  --axis_stroke_color '#374151' \
  --axis_stroke_width 4 \
  --line_stroke_color '#9CA3AF' \
  --line_stroke_width 1.5 \
  --species '<i>Homo sapiens</i>' \
  --plot_title 'Human mitochondrial feature presentation' \
  --plot_title_position top \
  --legend right \
  -o cli_feature_presentation \
  -f svg
```
<!-- executable:H-CLI-11:end -->

<!-- executable:H-CLI-01:start -->
```bash
gbdraw linear \
  --gbk NC_001416.gb \
  --show_labels none \
  --scale_style ruler \
  --track_layout above \
  --ruler_on_axis \
  -o lambda_genbank \
  -f svg

gbdraw linear \
  --gff NC_001416.gff3 \
  --fasta NC_001416.fna \
  --show_labels none \
  --scale_style ruler \
  --track_layout above \
  --ruler_on_axis \
  -o lambda_gff3 \
  -f svg
```
<!-- executable:H-CLI-01:end -->

<!-- executable:H-CLI-02:start -->
```bash
gbdraw circular \
  --records_table tables/records.tsv \
  --multi_record_canvas \
  --multi_record_size_mode equal \
  --circular_track_order ticks,features,gc_content,gc_skew \
  --circular_track_axis_index 1 \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels out \
  --label_font_size 10 \
  --plot_title 'Four complete metazoan mitochondrial genomes' \
  --plot_title_position bottom \
  --legend right \
  -o record_table \
  -f svg

gbdraw linear \
  --gbk NC_001416.gb NC_042057.1.gb \
  --record_id NC_001416.1 \
  --record_id NC_042057.1 \
  --comparisons_table tables/comparisons.tsv \
  --bitscore 50 \
  --evalue 1e-2 \
  --identity 0 \
  --alignment_length 0 \
  --pairwise_match_style curve \
  --show_labels none \
  --track_layout above \
  --scale_style ruler \
  --ruler_on_axis \
  --legend right \
  -o comparison_table \
  -f svg

gbdraw circular \
  --gbk HmmtDNA.gbk \
  --conservation_table tables/conservation.tsv \
  --conservation_reference subject \
  --bitscore 50 \
  --evalue 1e-5 \
  --identity 40 \
  --alignment_length 50 \
  --conservation_ring_width 18 \
  --conservation_ring_gap 4 \
  --track_type middle \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels out \
  --plot_title 'Complete metazoan mitochondrial TLOSATX evidence' \
  --plot_title_position top \
  --legend right \
  -o conservation_table \
  -f svg

gbdraw circular \
  --gbk NC_001879.gbk \
  --annotation_table nicotiana-tabacum-regions.tsv \
  -d modified_default_colors.tsv \
  --qualifier_priority qualifier_priority.tsv \
  --circular_track_slot 'ticks:ticks@side=outside,tick_label_layout=label_out_tick_in' \
  --circular_track_slot 'features:features@side=overlay,lane_direction=split,w=48px' \
  --circular_track_slot 'plastome_regions:annotations@set_id=plastome_regions,side=inside,w=30px,show_labels=true,padding_px=1,overflow=compress' \
  --circular_track_slot 'gc_content:dinucleotide_content@side=inside,w=36px,nt=GC,legend_label=GC content' \
  --circular_track_slot 'gc_skew:dinucleotide_skew@side=inside,w=32px,nt=GC,legend_label=GC skew' \
  --separate_strands \
  --track_type tuckin \
  --labels none \
  --plot_title '<i>Nicotiana tabacum</i> plastome regions' \
  --plot_title_position top \
  --legend right \
  -o annotation_table \
  -f svg

gbdraw circular \
  --gbk HmmtDNA.gbk \
  --circular_track_table tables/tracks.tsv \
  --track_type middle \
  --window 500 \
  --step 50 \
  --species '<i>Homo sapiens</i>' \
  --qualifier_priority cds_gene_qualifier_priority.tsv \
  --labels out \
  --legend left \
  -o track_table \
  -f svg
```
<!-- executable:H-CLI-02:end -->

<!-- executable:H-PY-03:start -->
```python
from pathlib import Path

from pandas import DataFrame

from gbdraw import (
    CircularLayout,
    CircularOptions,
    CircularTrackOptions,
    DepthTrackOptions,
    Diagram,
    FeatureOptions,
    LabelOptions,
    TitleOptions,
    draw_circular,
    read_genbank,
)
from gbdraw.api import AnnotationOptions, CircularTrackSlot, ScalarSpec


records = read_genbank(
    [Path("AP027133.gb"), Path("NC_001879.gbk"), Path("HmmtDNA.gbk")]
)
assert [(record.id, len(record)) for record in records] == [
    ("AP027133.1", 606_194),
    ("NC_001879.2", 155_943),
    ("NC_012920.1", 16_569),
]
assert all(record.annotations.get("topology") == "circular" for record in records)

label_whitelist = DataFrame(
    [
        [
            "CDS",
            "gene",
            "^(rpoB|secA|polC|rbcL|psaA|atpB|ndhF|"
            "ND1|COX1|ATP6|CYTB)$",
        ]
    ],
    columns=["feature_type", "qualifier", "keyword"],
)

track_slots = (
    CircularTrackSlot(
        id="features",
        renderer="features",
        side="overlay",
        params={"lane_direction": "split"},
    ),
    CircularTrackSlot(
        id="plastome_regions",
        renderer="annotations",
        side="inside",
        radius=ScalarSpec(0.70),
        width=ScalarSpec(24, "px"),
        params={
            "set_id": "plastome_regions",
            "show_labels": True,
            "padding_px": 1,
            "overflow": "compress",
        },
        inner_gap_px=2,
        outer_gap_px=2,
    ),
    CircularTrackSlot(
        id="depth",
        renderer="depth",
        side="inside",
        radius=ScalarSpec(0.59),
        width=ScalarSpec(0.10),
        params={"track_index": 0, "legend_label": "DRR394922 mean depth"},
        inner_gap_px=3,
        outer_gap_px=3,
    ),
    CircularTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        side="inside",
        radius=ScalarSpec(0.45),
        width=ScalarSpec(0.08),
        params={"nt": "GC", "legend_label": "GC content"},
        inner_gap_px=3,
        outer_gap_px=3,
    ),
)

options = CircularOptions(
    features=FeatureOptions(
        types=("CDS", "rRNA", "tRNA", "tmRNA", "repeat_region"),
        default_colors=Path("modified_default_colors.tsv"),
    ),
    labels=LabelOptions(
        whitelist=label_whitelist,
        qualifier_priority=Path("qualifier_priority.tsv"),
    ),
    annotations=AnnotationOptions(table_file="nicotiana-tabacum-regions.tsv"),
    tracks=CircularTrackOptions(slots=track_slots, center_reserved_radius=56),
    depth_tracks=(
        DepthTrackOptions(
            source=(Path("AP027133.DRR394922.depth-1kb.tsv"), None, None),
            label="DRR394922 mean depth",
            color="#2563eb",
            large_tick_interval=20,
            small_tick_interval=10,
        ),
    ),
    depth_window=1,
    depth_step=1000,
    title=TitleOptions(
        text="Genome tracks and annotations",
        position="top",
        font_size=26,
    ),
    legend="right",
    keep_full_definition_with_title=True,
    config_overrides={
        "canvas.strandedness": True,
        "canvas.resolve_overlaps": True,
        "canvas.circular.track_type": "tuckin",
        "labels.circular.scope": "outer",
        "labels.circular.placement": "horizontal",
        "labels.font_size.long": 11,
        "objects.features.block_stroke_color": "#4b5563",
        "objects.features.block_stroke_width.long": 0.8,
        "objects.features.line_stroke_color": "#9ca3af",
        "objects.features.line_stroke_width.long": 1.5,
    },
)

tracks_diagram = draw_circular(
    records,
    options=options,
    layout=CircularLayout(size="equal", positions=("#1@1", "#2@1", "#3@2")),
)
tracks_svg = tracks_diagram.to_svg()
tracks_bytes = tracks_diagram.to_bytes("svg")
tracks_path = tracks_diagram.save(Path("python_tracks_annotations.svg"))

assert isinstance(tracks_diagram, Diagram)
assert tracks_diagram.mode == "circular"
assert tracks_svg.encode("utf-8") == tracks_bytes == tracks_path.read_bytes()
print("Saved python_tracks_annotations.svg")
```
<!-- executable:H-PY-03:end -->

<!-- executable:H-PY-05:start -->
```python
from datetime import datetime, timezone
from pathlib import Path
from tempfile import TemporaryDirectory

from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    CircularOutputOptions,
    CircularRequestPlan,
    CircularRequestTrackOptions,
    InMemoryRecordSource,
    PreparedDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    RequestRenderResult,
    SessionConversionError,
    SessionDocument,
    build_request_plan_diagram,
    load_gbks,
    load_session_document,
    materialize_session,
    plan_request,
    render_prepared_request,
    render_session,
    save_session_document,
    session_to_request,
)


record = load_gbks(["HmmtDNA.gbk"])[0]
request = CircularDiagramRequest(
    records=(
        RecordInput(
            source=InMemoryRecordSource(record),
            record_key="human-mitochondrion",
        ),
    ),
    options=CircularDiagramOptions(
        tracks=CircularRequestTrackOptions(),
        output=CircularOutputOptions(legend="right"),
    ),
    output=RenderOutputRequest(
        output_prefix="typed_request",
        formats=("svg",),
    ),
)

plan = plan_request(request)
assert isinstance(plan, CircularRequestPlan)
assert plan.mode == "circular"
plan.preflight_outputs()

prepared = build_request_plan_diagram(plan)
assert isinstance(prepared, PreparedDiagramRequest)

result = render_prepared_request(prepared)
assert isinstance(result, RequestRenderResult)
assert result.mode == "circular"
assert tuple(path.name for path in result.output_paths) == ("typed_request.svg",)
rendered_svg_bytes = result.output_paths[0].read_bytes()

session_path = Path("typed_request.session.json")
session_document = save_session_document(
    session_path,
    request,
    title="Typed human mitochondrial diagram",
    created_at=datetime(2026, 8, 3, tzinfo=timezone.utc),
)
assert isinstance(session_document, SessionDocument)

loaded_document = load_session_document(session_path)
assert loaded_document.mode == "circular"
assert loaded_document.has_canonical_request

with TemporaryDirectory(prefix="gbdraw-session-replay-") as replay_directory:
    with materialize_session(
        loaded_document,
        output_directory=replay_directory,
    ) as materialized:
        replay_request = session_to_request(materialized)
        assert isinstance(replay_request, CircularDiagramRequest)
        replay_result = render_session(materialized)
        assert isinstance(replay_result, RequestRenderResult)
        replayed_svg_bytes = replay_result.output_paths[0].read_bytes()
    resources_expired = not materialized.active
replay_output_removed = not Path(replay_directory).exists()

assert rendered_svg_bytes == replayed_svg_bytes
assert resources_expired
assert replay_output_removed

wrong_mode_payload = loaded_document.to_dict()
wrong_mode_payload["renderRequest"]["diagramOptions"]["tracks"][
    "linearTrackAxisIndex"
] = 0
wrong_mode_document = load_session_document(wrong_mode_payload)
wrong_mode_error = ""

with TemporaryDirectory(prefix="gbdraw-wrong-mode-") as invalid_output_directory:
    with materialize_session(
        wrong_mode_document,
        output_directory=invalid_output_directory,
    ) as invalid_materialized:
        try:
            session_to_request(invalid_materialized)
        except SessionConversionError as error:
            wrong_mode_error = str(error)
        else:
            raise AssertionError("The wrong-mode session was accepted")

assert "Circular request cannot contain Linear track values" in wrong_mode_error

print("Rendered typed_request.svg")
print("Replayed the current Circular session")
print(f"Rejected wrong-mode session: {wrong_mode_error}")
```
<!-- executable:H-PY-05:end -->

<!-- executable:H-PY-01:start -->
```python
from pathlib import Path

from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw import (
    CircularOptions,
    Diagram,
    LabelOptions,
    TitleOptions,
    draw_circular,
    read_genbank,
)


cds_label_whitelist = DataFrame(
    [("CDS", "gene", ".+")],
    columns=["feature_type", "qualifier", "keyword"],
)
labels = LabelOptions(
    whitelist=cds_label_whitelist,
    qualifier_priority=Path("cds_gene_qualifier_priority.tsv"),
)
label_overrides = {
    "labels.circular.scope": "outer",
    "labels.circular.placement": "horizontal",
    "labels.font_size.short": 18,
    "labels.font_size.long": 18,
}


circular_record = read_genbank(Path("HmmtDNA.gbk"))[0]
assert isinstance(circular_record, SeqRecord)
assert (circular_record.id, len(circular_record)) == ("NC_012920.1", 16_569)

circular_options = CircularOptions(
    labels=labels,
    config_overrides=label_overrides,
)
circular_diagram = draw_circular(circular_record, options=circular_options)
circular_svg = circular_diagram.to_svg()
circular_bytes = circular_diagram.to_bytes("svg")
circular_path = circular_diagram.save(Path("python_circular.svg"))

assert isinstance(circular_diagram, Diagram)
assert circular_diagram.mode == "circular"
assert circular_svg.encode("utf-8") == circular_bytes
assert circular_path.read_bytes() == circular_bytes

multi_records = read_genbank(
    [
        Path("HmmtDNA.gbk"),
        Path("NC_002333.2.gb"),
        Path("NC_024511.2.gb"),
        Path("NC_001328.1.gb"),
    ]
)
assert [(record.id, len(record)) for record in multi_records] == [
    ("NC_012920.1", 16_569),
    ("NC_002333.2", 16_596),
    ("NC_024511.2", 19_524),
    ("NC_001328.1", 13_794),
]
assert all(record.annotations.get("topology") == "circular" for record in multi_records)

multi_record_options = CircularOptions(
    labels=labels,
    title=TitleOptions(
        text="Complete metazoan mitochondrial genomes",
        position="top",
    ),
    keep_full_definition_with_title=True,
    config_overrides=label_overrides,
)
multi_record_diagram = draw_circular(multi_records, options=multi_record_options)
multi_record_svg = multi_record_diagram.to_svg()
multi_record_bytes = multi_record_diagram.to_bytes("svg")
multi_record_path = multi_record_diagram.save(Path("python_multi_record.svg"))

assert isinstance(multi_record_diagram, Diagram)
assert multi_record_diagram.mode == "circular"
assert len(multi_record_diagram.records) == 4
assert multi_record_svg.encode("utf-8") == multi_record_bytes
assert multi_record_path.read_bytes() == multi_record_bytes

print("Saved python_circular.svg and python_multi_record.svg")
```
<!-- executable:H-PY-01:end -->

<!-- executable:H-PY-02:start -->
```python
from pathlib import Path

from Bio.SeqRecord import SeqRecord

from gbdraw import (
    Diagram,
    FeatureOptions,
    LinearComparisonOptions,
    LinearOptions,
    Thresholds,
    TitleOptions,
    draw_linear,
    read_genbank,
)


comparison_records = read_genbank(
    [Path("NC_001416.gb"), Path("NC_042057.1.gb")]
)
assert [(record.id, len(record)) for record in comparison_records] == [
    ("NC_001416.1", 48_502),
    ("NC_042057.1", 42_925),
]
assert all(isinstance(record, SeqRecord) for record in comparison_records)

comparison_options = LinearComparisonOptions(
    blast_files=("lambda-de3.losatn.tsv",),
)
comparison_diagram = draw_linear(
    comparison_records,
    options=LinearOptions(
        features=FeatureOptions(types=("CDS",)),
        comparisons=comparison_options,
        thresholds=Thresholds(evalue=1e-5),
        title=TitleOptions(text="Whole-genome Lambda and DE3 LOSATN comparison", position="top"),
        legend="none",
    ),
)
comparison_svg = comparison_diagram.to_svg()
comparison_bytes = comparison_diagram.to_bytes("svg")
comparison_path = comparison_diagram.save(Path("python_linear_comparison.svg"))

assert isinstance(comparison_diagram, Diagram)
assert comparison_diagram.mode == "linear"
assert comparison_svg.encode("utf-8") == comparison_bytes
assert comparison_path.read_bytes() == comparison_bytes

print("Saved python_linear_comparison.svg")
```
<!-- executable:H-PY-02:end -->

<!-- executable:H-PY-04:start -->
```python
from io import StringIO
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gbdraw import Diagram, draw_circular, draw_linear, read_gff


gff_records = read_gff(
    Path("NC_001416.gff3"),
    Path("NC_001416.fna"),
    features=("CDS", "gene"),
)
assert len(gff_records) == 1
gff_record = gff_records[0]
assert isinstance(gff_record, SeqRecord)
assert (gff_record.id, len(gff_record)) == ("NC_001416.1", 48_502)

gff_diagram = draw_linear(gff_record)
gff_svg = gff_diagram.to_svg()
gff_bytes = gff_diagram.to_bytes("svg")
gff_path = Path("python_gff3.svg")
gff_path.write_bytes(gff_bytes)

assert isinstance(gff_diagram, Diagram)
assert gff_diagram.mode == "linear"
assert gff_svg.encode("utf-8") == gff_bytes
assert gff_path.read_bytes() == gff_bytes

genbank_text = Path("HmmtDNA.gbk").read_text(encoding="utf-8")
memory_record = SeqIO.read(StringIO(genbank_text), "genbank")
assert isinstance(memory_record, SeqRecord)
assert (memory_record.id, len(memory_record)) == ("NC_012920.1", 16_569)

memory_diagram = draw_circular(memory_record)
memory_svg = memory_diagram.to_svg()
memory_bytes = memory_diagram.to_bytes("svg")
memory_path = Path("python_memory.svg")
memory_path.write_bytes(memory_bytes)

assert isinstance(memory_diagram, Diagram)
assert memory_diagram.mode == "circular"
assert memory_svg.encode("utf-8") == memory_bytes
assert memory_path.read_bytes() == memory_bytes

print("Wrote python_gff3.svg and python_memory.svg")
```
<!-- executable:H-PY-04:end -->
