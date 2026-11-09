"""Feature queries and feature density tracks: matching, window values, drawing, Legend rows, request round trip."""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.feature_density import (
    FeatureDensityTrack,
    feature_density_layers,
    feature_density_values,
    parse_feature_density_words,
    read_feature_density_table,
)
from gbdraw.api.assembly_inputs import CircularAssemblyInputs, LinearAssemblyInputs
from gbdraw.api.diagram import assemble_circular_diagram_from_record, assemble_linear_diagram_from_records
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.features.query import FeatureQuery


def _feature(start: int, end: int, feature_type: str = "gene", **qualifiers: str) -> SeqFeature:
    return SeqFeature(
        SimpleLocation(start, end, strand=1),
        type=feature_type,
        qualifiers={key: [value] for key, value in qualifiers.items()},
    )


def _join(parts: list[tuple[int, int]], feature_type: str = "gene") -> SeqFeature:
    return SeqFeature(CompoundLocation([SimpleLocation(s, e, strand=1) for s, e in parts]), type=feature_type)


def _record(features: list[SeqFeature], length: int = 100, topology: str = "circular") -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * (length // 4)), id="r1", name="r1", description="test record")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = topology
    record.features = [SeqFeature(SimpleLocation(0, length, strand=1), type="source"), *features]
    return record


# A 100 bp record: genes at 5-15, 15-25 (touching), 20-50 (overlapping), a join 92-98 + 0-4 across the origin,
# and a join 60-65 + 68-70 whose two parts fall in one window.
_FEATURES = [
    _feature(5, 15, gene="dnaA"),
    _feature(15, 25, gene="dnaN"),
    _feature(20, 50, gene="recF", note="long gene"),
    _join([(92, 100), (0, 4)]),
    _join([(60, 65), (68, 70)]),
]


def _values(metric: str, *, wraps: bool, window: int = 10, step: int = 10) -> list[float]:
    track = FeatureDensityTrack("t", FeatureQuery("gene"), metric=metric)  # type: ignore[arg-type]
    return feature_density_values(_record(_FEATURES), track, window=window, step=step, wraps=wraps).values.tolist()


# --- FeatureQuery -------------------------------------------------------------------------------------------------


def test_query_matches_exact_type_and_star() -> None:
    gene, cds = _feature(0, 10, "gene"), _feature(0, 10, "CDS")
    assert FeatureQuery("gene").matches(gene)
    assert not FeatureQuery("gene").matches(cds)
    assert FeatureQuery("*").matches(gene) and FeatureQuery("*").matches(cds)


def test_query_searches_qualifier_values_case_insensitively() -> None:
    ltr = _feature(0, 10, "repeat_region", rpt_type="long_terminal_repeat")
    assert FeatureQuery("repeat_region", "rpt_type", "LONG_TERMINAL").matches(ltr)
    assert FeatureQuery("repeat_region", "RPT_TYPE", "^long").matches(ltr)
    assert not FeatureQuery("repeat_region", "rpt_type", "^terminal").matches(ltr)
    # A qualifier without a value matches the features that carry it.
    assert FeatureQuery("*", "rpt_type").matches(ltr)
    assert not FeatureQuery("*", "gene").matches(ltr)


def test_query_reads_pseudo_qualifiers() -> None:
    feature = _feature(4, 10, "gene")
    assert FeatureQuery("gene", "location", r"^4\.\.10$").matches(feature, record_id="r1")
    assert FeatureQuery("gene", "record_location", "^r1:").matches(feature, record_id="r1")


def test_query_length_is_the_bases_of_every_part() -> None:
    joined = _join([(0, 10), (50, 60)])  # spans 60 bp, covers 20
    assert FeatureQuery("gene", min_length=20).matches(joined)
    assert not FeatureQuery("gene", min_length=21).matches(joined)
    assert FeatureQuery("gene", max_length=20).matches(joined)
    assert not FeatureQuery("gene", min_length=5, max_length=19).matches(joined)


def test_query_refuses_unusable_conditions() -> None:
    with pytest.raises(ParseError, match="Invalid regex"):
        FeatureQuery("gene", "gene", "(")
    with pytest.raises(ValidationError, match="needs a qualifier"):
        FeatureQuery("gene", None, "dnaA")
    with pytest.raises(ValidationError, match="above its max_length"):
        FeatureQuery("gene", min_length=10, max_length=5)
    with pytest.raises(ValidationError):
        FeatureQuery("gene", min_length=-1)
    with pytest.raises(ValidationError):
        FeatureQuery("")


def test_query_describes_itself_in_words() -> None:
    assert FeatureQuery("*").describe() == "any feature"
    assert (
        FeatureQuery("repeat_region", "rpt_type", "long_terminal_repeat").describe()
        == "repeat_region · rpt_type = long_terminal_repeat"
    )
    assert FeatureQuery("gene", min_length=5000).describe() == "gene · length ≥ 5 kbp"
    assert FeatureQuery("gene", min_length=500, max_length=2_500_000).describe() == "gene · length 500 bp–2.5 Mbp"


# --- Window values ------------------------------------------------------------------------------------------------


def test_count_on_a_circular_record_wraps_round_the_origin() -> None:
    # Windows centre on 0, 10, ..., 90 and span centre-5 .. centre+5; the first reads 95..100 and 0..5.
    assert _values("count", wraps=True) == [1, 1, 2, 1, 1, 1, 1, 1, 0, 1]


def test_count_on_a_linear_record_stops_at_its_ends() -> None:
    # The first window reads 0..5 only; the join across the origin is in it by its part 0..4.
    assert _values("count", wraps=False) == [1, 1, 2, 1, 1, 1, 1, 1, 0, 1]
    # A feature at 96..100 is in the first window only when it wraps.
    track = FeatureDensityTrack("t", FeatureQuery("gene"))
    record = _record([_feature(96, 100)])
    wrapped = feature_density_values(record, track, window=10, step=10, wraps=True).values.tolist()
    clipped = feature_density_values(record, track, window=10, step=10, wraps=False).values.tolist()
    assert wrapped[0] == 1 and clipped[0] == 0
    assert wrapped[1:] == clipped[1:] == [0] * 9


def test_a_join_counts_once_and_a_long_feature_counts_in_every_window_it_crosses() -> None:
    track = FeatureDensityTrack("t", FeatureQuery("gene"))
    record = _record([_join([(60, 65), (68, 70)]), _feature(10, 90)])
    values = feature_density_values(record, track, window=10, step=10, wraps=False).values.tolist()
    # Window 55..65 and 65..75 hold one part each of the join; every window from 5..15 to 85..95 holds the long gene.
    assert values == [0, 1, 1, 1, 1, 1, 2, 2, 1, 1]


def test_coverage_is_the_percentage_of_the_window_in_the_union_of_parts() -> None:
    # The first window reads 95..100 and 0..5: 9 of its 10 bases are in the join across the origin;
    # touching and overlapping genes count each base once (15..25 and 20..50).
    assert _values("coverage", wraps=True) == [90.0, 100.0, 100.0, 100.0, 100.0, 50.0, 50.0, 20.0, 0.0, 30.0]
    # Clipped at the ends: the first window reads 0..5, 4 of its 5 bases covered.
    assert _values("coverage", wraps=False) == [80.0, 100.0, 100.0, 100.0, 100.0, 50.0, 50.0, 20.0, 0.0, 30.0]


def _reference_values(record: SeqRecord, metric: str, *, wraps: bool, window: int, step: int) -> list[float]:
    """The window values base by base: the windows the GC track reads, each feature as the set of its bases."""

    from gbdraw.analysis.windows import sliding_windows

    length = len(record.seq)
    window = min(window, length) if wraps else window
    windows = sliding_windows(record, window, step, circular=wraps)
    features = [
        {base for part in feature.location.parts for base in range(int(part.start), int(part.end))}
        for feature in record.features
        if feature.type != "source"
    ]
    covered = set().union(*features)
    values = []
    for start, end in zip(windows.starts.tolist(), windows.ends.tolist()):
        if wraps:
            bases = {base % length for base in range(start, end)}
        else:
            start, end = max(0, start), min(length, end)
            bases = set(range(start, end))
        if metric == "count":
            values.append(float(sum(1 for feature in features if feature & bases)))
        else:
            values.append(100.0 * len(bases & covered) / max(end - start, 1))
    return values


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("wraps", [True, False])
@pytest.mark.parametrize("metric", ["count", "coverage"])
def test_window_values_equal_a_base_by_base_reference(seed: int, wraps: bool, metric: str) -> None:
    rng = np.random.default_rng(seed)
    length = 240
    features = []
    for _ in range(25):
        parts = sorted({(int(start), int(min(length, start + rng.integers(1, 40)))) for start in rng.integers(0, length, rng.integers(1, 4))})
        features.append(_join(parts) if len(parts) > 1 else _feature(*parts[0]))
    # A feature across the origin and one of the whole record.
    features += [_join([(length - 7, length), (0, 5)]), _feature(0, length)]
    record = _record(features, length=length)
    track = FeatureDensityTrack("t", FeatureQuery("gene"), metric=metric)  # type: ignore[arg-type]
    for window, step in ((20, 10), (35, 15), (240, 60), (500, 100)):
        values = feature_density_values(record, track, window=window, step=step, wraps=wraps).values.tolist()
        assert values == pytest.approx(_reference_values(record, metric, wraps=wraps, window=window, step=step))


def test_a_track_window_and_step_override_the_default() -> None:
    track = FeatureDensityTrack("t", FeatureQuery("gene"), window=50, step=25)
    values = feature_density_values(_record(_FEATURES), track, window=10, step=10, wraps=True)
    assert values.samples.tolist() == [0, 25, 50, 75]
    # A window may not read more than the whole circle.
    huge = FeatureDensityTrack("t", FeatureQuery("gene"), window=500, step=50)
    assert feature_density_values(_record(_FEATURES), huge, window=10, step=10, wraps=True).values.tolist() == [5, 5]


def test_layers_share_one_scale_across_records() -> None:
    track = FeatureDensityTrack("t", FeatureQuery("gene"))
    dense = _record([_feature(i, i + 2) for i in range(0, 100, 3)])
    sparse = _record([_feature(0, 2)])
    layers = feature_density_layers([dense, sparse], [track], windows=[(10, 10), (10, 10)], wraps=[True, True])
    assert layers[0]["t"].scale_max == layers[1]["t"].scale_max == max(layers[0]["t"].values.values)


def test_source_features_are_never_counted() -> None:
    track = FeatureDensityTrack("t", FeatureQuery("*"))
    values = feature_density_values(_record([]), track, window=10, step=10, wraps=True)
    assert np.all(values.values == 0)


# --- Track definitions --------------------------------------------------------------------------------------------


def test_track_words_and_table_define_the_same_track(tmp_path: Path) -> None:
    words = parse_feature_density_words(
        ["id=ltr", "type=repeat_region", "qualifier=rpt_type", "value=long_terminal_repeat", "metric=coverage",
         "style=histogram", "window=20000", "step=5000", "color=#d08a3a", "caption=LTR"],
        where="test",
    )
    table = tmp_path / "density.tsv"
    table.write_text(
        "id\ttype\tqualifier\tvalue\tmetric\tstyle\twindow\tstep\tcolor\tcaption\n"
        "ltr\trepeat_region\trpt_type\tlong_terminal_repeat\tcoverage\thistogram\t20000\t5000\t#d08a3a\tLTR\n"
    )
    assert read_feature_density_table(str(table)) == (words,)


@pytest.mark.parametrize(
    ("words", "message"),
    [
        (["type=gene"], "needs an id"),
        (["id=a"], "needs a feature type"),
        (["id=a", "type=gene", "metric=sum"], "metric must be one of"),
        (["id=a", "type=gene", "window=0"], "window must be a whole number"),
        (["id=a b", "type=gene"], "track id"),
        (["id=a", "type=gene", "shape=x"], "unknown feature density field"),
        (["id=a", "type=gene", "color=notacolor"], "notacolor"),
    ],
)
def test_bad_track_words_are_refused(words: list[str], message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        parse_feature_density_words(words, where="test")


def test_two_tracks_of_one_id_are_refused() -> None:
    track = FeatureDensityTrack("a", FeatureQuery("gene"))
    with pytest.raises(ValidationError, match="Two feature density tracks"):
        CircularDiagramOptions(feature_density_tracks=(track, track))


# --- Drawing ------------------------------------------------------------------------------------------------------


@pytest.fixture(scope="module")
def cfg() -> GbdrawConfig:
    return GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))


_TRACKS = (
    FeatureDensityTrack("genes", FeatureQuery("gene")),
    FeatureDensityTrack("long", FeatureQuery("gene", min_length=20), metric="coverage", style="histogram",
                        color="#d08a3a", caption="Long genes"),
)


def _slot_group(svg: str, slot_id: str) -> str:
    match = re.search(rf'<g[^>]*(?:data-gbdraw-slot-id|id)="{slot_id}"[^>]*>(.*?)</g>', svg, re.S)
    assert match is not None, f"no group for {slot_id}"
    return match.group(1)


def test_circular_draws_a_heatmap_band_and_a_histogram_with_legend_rows(cfg: GbdrawConfig) -> None:
    record = _record(_FEATURES, length=1000)
    svg = assemble_circular_diagram_from_record(
        record, CircularAssemblyInputs(cfg=cfg, feature_density_tracks=_TRACKS)
    ).tostring()
    heatmap = _slot_group(svg, "density_genes")
    assert heatmap.count("<path") >= 1
    assert 'data-gbdraw-slot-renderer="feature_density"' in svg
    assert "density_long" in svg
    # The heatmap's Legend row is a gradient from 0 to the scale's maximum; the histogram's names its maximum.
    assert "gene per" in svg
    assert "Long genes" in svg


def test_linear_draws_a_heatmap_band_and_a_histogram(cfg: GbdrawConfig) -> None:
    record = _record(_FEATURES, length=1000, topology="linear")
    drawing = assemble_linear_diagram_from_records([record], LinearAssemblyInputs(cfg=cfg, feature_density_tracks=_TRACKS))
    svg = getattr(drawing, "drawing", drawing).tostring()
    assert _slot_group(svg, "density_genes").count("<path") + _slot_group(svg, "density_genes").count("<rect") >= 1
    assert "density_long" in svg
    assert "gene per" in svg and "Long genes" in svg


def _slot_y_extent(svg: str, slot_id: str) -> tuple[float, float]:
    """The top and bottom of a Linear slot's path ink on the canvas: the slot's own vertical translate plus the
    y of every path point."""

    match = re.search(rf'<g[^>]*data-gbdraw-slot-id="{slot_id}"[^>]*transform="([^"]*)"[^>]*>(.*?)</g>', svg, re.S)
    assert match is not None, f"no group for {slot_id}"
    origin = sum(float(y) for y in re.findall(r"translate\([^,)]+,\s*([-\d.]+)\)", match.group(1)))
    ys = [float(y) for y in re.findall(r"[ML ,]\s*-?[\d.]+[ ,](-?[\d.]+)", " ".join(re.findall(r' d="([^"]*)"', match.group(2))))]
    assert ys, f"no path in {slot_id}"
    return origin + min(ys), origin + max(ys)


def test_a_linear_histogram_keeps_its_height_in_the_slot_stack(cfg: GbdrawConfig) -> None:
    # The GC content slot below a histogram starts below it, as below a Depth slot.
    record = _record(_FEATURES, length=1000, topology="linear")
    # GC that swings between AT-rich and GC-rich stretches, so the GC content track reaches its full height.
    record.seq = Seq("".join(("AT" if (position // 100) % 2 else "GC") * 25 for position in range(0, 1000, 50)))
    drawing = assemble_linear_diagram_from_records(
        [record],
        LinearAssemblyInputs(
            cfg=cfg,
            feature_density_tracks=_TRACKS,
            linear_track_slots=["features:features", "long:feature_density@track=long", "gc:dinucleotide_content"],
        ),
    )
    svg = getattr(drawing, "drawing", drawing).tostring()
    histogram_top, histogram_bottom = _slot_y_extent(svg, "long")
    gc_top, _ = _slot_y_extent(svg, "gc")
    assert histogram_bottom > histogram_top
    assert gc_top >= histogram_bottom


def test_an_explicit_slot_names_its_track_and_an_unknown_track_is_refused(cfg: GbdrawConfig) -> None:
    record = _record(_FEATURES, length=1000)
    svg = assemble_circular_diagram_from_record(
        record,
        CircularAssemblyInputs(
            cfg=cfg,
            feature_density_tracks=_TRACKS,
            circular_track_slots=["features:features", "g:feature_density@track=genes", "ticks:ticks"],
        ),
    ).tostring()
    assert 'data-gbdraw-slot-id="g"' in svg
    assert "density_long" not in svg and "Long genes" not in svg
    with pytest.raises(ValidationError, match="nope"):
        assemble_circular_diagram_from_record(
            record,
            CircularAssemblyInputs(
                cfg=cfg, feature_density_tracks=_TRACKS, circular_track_slots=["d:feature_density@track=nope"]
            ),
        )


# --- Request ------------------------------------------------------------------------------------------------------


DENSITY_SESSION = Path(__file__).resolve().parent / "fixtures/sessions/HmmtDNA-feature-density.v51-cli.gbdraw-session.json.gz"


def test_a_cli_session_with_feature_density_tracks_replays_as_it_was_drawn(tmp_path: Path) -> None:
    from gbdraw.api import load_session_document
    from gbdraw.circular import circular_main
    from gbdraw.session import session_drawing_artifacts

    drawing = session_drawing_artifacts(load_session_document(DENSITY_SESSION), "circular").fields
    stored = drawing["renderRequest"]["diagramOptions"]["featureDensityTracks"]
    assert [track["id"] for track in stored] == ["genes", "trna"]

    circular_main(["--session", str(DENSITY_SESSION), "-o", str(tmp_path / "replay"), "-f", "svg"])

    svg = (tmp_path / "replay.svg").read_text()
    assert 'data-gbdraw-slot-id="density_genes"' in svg and 'data-gbdraw-slot-id="density_trna"' in svg
    assert "gene per 1 kbp" in svg and "tRNA, % of 1 kbp" in svg
