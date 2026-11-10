"""A crop's windows read the bases beside it from the record it was cut from (OV-2890); a Session keeps them.

A crop given in memory, or a crop of a file that is no longer at hand, is written to a Session as the record it was
cut from and the crop's window in it, so the replay crops it again and draws what the live figure drew (D-64).
"""

from __future__ import annotations

import copy
import random
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    CircularBatchRequest,
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    apply_region_specs,
    load_session_document,
    materialize_session,
    parse_region_specs,
    render_request,
    render_session,
    resolve_request,
    save_session_document,
)
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions, LinearMultiRecordOptions
from gbdraw.api.request_render import _materialized_record_inputs, resolve_request_records
from gbdraw.api.requests import GenBankInputSource, RecordCollectionOptions, RecordPresentation, RecordSelector
from gbdraw.session_container import read_session_container
from gbdraw.io.regions import RecordWindow
from gbdraw.core.record_metadata import held_record_origin

_TRACKS = {"canvas.show_gc": True, "canvas.show_skew": True}


def _sequence(seed: int) -> str:
    """GC-rich flanks round an AT-rich middle: a window that reads a flank shows it."""

    rng = random.Random(seed)
    return "".join(rng.choice(bases) for bases in ("GC", "AT", "GC") for _ in range(20_000))


def _record(record_id: str, topology: str, seed: int = 3) -> SeqRecord:
    record = SeqRecord(
        Seq(_sequence(seed)), id=record_id, name=record_id, annotations={"topology": topology, "molecule_type": "DNA"}
    )
    record.features = [SeqFeature(FeatureLocation(25_000, 26_000, 1), type="CDS", qualifiers={"product": ["x"]})]
    return record


def _replayed(request, tmp_path: Path) -> tuple[str, str]:
    """The live figure and the figure of its saved Session."""

    render_request(request)
    live = (tmp_path / "live.svg").read_text()
    path = tmp_path / "saved.gbdraw"
    save_session_document(path, request)
    replay_dir = tmp_path / "replay"
    replay_dir.mkdir()
    with materialize_session(load_session_document(path), output_directory=replay_dir) as materialized:
        render_session(materialized)
    return live, next(replay_dir.glob("*.svg")).read_text()


def _output(tmp_path: Path) -> RenderOutputRequest:
    return RenderOutputRequest(output_prefix="live", output_directory=tmp_path, overwrite=True)


@pytest.mark.parametrize(
    "mode,topology,spec,reverse",
    [
        ("linear", "linear", "src:20001-40000", False),
        ("linear", "linear", "src:20001-40000:rc", False),
        ("linear", "linear", "src:20001-40000", True),
        ("linear", "linear", "src:20001-40000:rc", True),
        ("circular", "circular", "src:50001-10000", False),
        ("circular", "linear", "src:20001-40000", False),
    ],
    ids=["linear", "linear-rc", "linear-reversed", "linear-rc-reversed", "circular-across-origin", "circular-linear"],
)
def test_a_session_redraws_an_in_memory_crop_with_its_flanks(
    mode: str, topology: str, spec: str, reverse: bool, tmp_path: Path
) -> None:
    crop = apply_region_specs([_record("src", topology)], parse_region_specs([spec]))[0]
    record = RecordInput(
        source=InMemoryRecordSource(crop), record_key="r", presentation=RecordPresentation(reverse_complement=reverse)
    )
    request = (
        LinearDiagramRequest(records=(record,), options=LinearDiagramOptions(config_overrides=_TRACKS),
                             output=_output(tmp_path))
        if mode == "linear"
        else CircularDiagramRequest(records=(record,), options=CircularDiagramOptions(config_overrides=_TRACKS),
                                    output=_output(tmp_path))
    )
    live, replay = _replayed(request, tmp_path)
    assert replay == live


def _saved_records(tmp_path: Path) -> list[dict]:
    document = read_session_container(tmp_path / "saved.gbdraw").document
    return document["drawings"][0]["renderRequest"]["records"]


def _circular_grid_or_batch(grouping: str, records: tuple[RecordInput, ...], tmp_path: Path):
    options = CircularDiagramOptions(config_overrides=_TRACKS)
    if grouping == "grid":
        return CircularDiagramRequest(records=records, options=options, output=_output(tmp_path), grouping="grid")
    return CircularBatchRequest(
        records=records, options=options,
        outputs=tuple(RenderOutputRequest(output_prefix=f"live-{index}", output_directory=tmp_path, overwrite=True)
                      for index in range(len(records))),
    )


@pytest.mark.parametrize("grouping", ["grid", "batch"])
def test_a_circular_grid_or_batch_session_redraws_in_memory_crops_with_their_flanks(
    grouping: str, tmp_path: Path
) -> None:
    """Crops given in memory to a Circular grid or batch, reversed after the cut or not, are saved as the records
    they were cut from and their ranges (PD-OI-044 scenario revision 5); a record the caller sliced itself has no
    source record and is saved as drawn."""

    crop = apply_region_specs([_record("src", "circular")], parse_region_specs(["src:20001-40000"]))[0]
    reversed_crop = apply_region_specs([_record("rev", "circular", 4)], parse_region_specs(["rev:15001-35000:rc"]))[0]
    sliced = _record("own", "circular", 5)[20_000:40_000]
    sliced.annotations = {"topology": "circular", "molecule_type": "DNA"}
    records = tuple(
        RecordInput(source=InMemoryRecordSource(record), record_key=key)
        for key, record in (("a", crop), ("b", reversed_crop), ("c", sliced))
    )
    request = _circular_grid_or_batch(grouping, records, tmp_path)
    render_request(request)
    path = tmp_path / "saved.gbdraw"
    save_session_document(path, request)
    saved = _saved_records(tmp_path)
    assert [(record["recordKey"], record["windows"]) for record in saved] == [
        ("a", [{"start": 20001, "end": 40000, "reverseComplement": False, "column": None}]),
        ("b", [{"start": 15001, "end": 35000, "reverseComplement": True, "column": None}]),
        ("c", []),
    ]
    replay_dir = tmp_path / "replay"
    replay_dir.mkdir()
    with materialize_session(load_session_document(path), output_directory=replay_dir) as materialized:
        render_session(materialized)
    live = sorted(tmp_path.glob("live*.svg"))
    assert len(live) == (1 if grouping == "grid" else 3)
    for drawn in live:
        assert (replay_dir / drawn.name).read_text() == drawn.read_text()


def test_two_crops_of_one_in_memory_record_store_one_copy_of_it(tmp_path: Path) -> None:
    """Two crops of one in-memory record (not windows of one record) read one resource of that record."""

    source = _record("src", "linear")
    first, second = apply_region_specs([source, source], parse_region_specs(["#1:5001-15000", "#2:30001-45000:rc"]))
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(first), record_key="a"),
                 RecordInput(InMemoryRecordSource(second), record_key="b")),
        options=LinearDiagramOptions(config_overrides=_TRACKS), layout=LinearMultiRecordOptions(),
        output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    saved = _saved_records(tmp_path)
    assert [record["windows"] for record in saved] == [
        [{"start": 5001, "end": 15000, "reverseComplement": False, "column": None}],
        [{"start": 30001, "end": 45000, "reverseComplement": True, "column": None}],
    ]
    assert saved[0]["source"] == saved[1]["source"]
    assert replay == live


def test_a_session_redraws_crops_of_a_file_no_longer_at_hand(tmp_path: Path) -> None:
    """A crop of each of two records of one file, given as the planner's records after the file is gone: each is
    written as the record it was cut from and its window."""

    path = tmp_path / "two.gb"
    SeqIO.write([_record("c1", "linear", 5), _record("c2", "linear", 6)], path, "genbank")
    options = LinearDiagramOptions(config_overrides=_TRACKS)
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(str(path)),
                selector=RecordSelector(raw="c1", record_id="c1", record_index=None),
                windows=(RecordWindow(35001, 45000, True),),
                record_key="a",
                presentation=RecordPresentation(grid_row=1, grid_column=1),
            ),
            RecordInput(
                source=GenBankInputSource(str(path)),
                selector=RecordSelector(raw="c2", record_id="c2", record_index=None),
                windows=(RecordWindow(15001, 40000, False),),
                record_key="b",
                presentation=RecordPresentation(grid_row=2, grid_column=1),
            ),
        ),
        options=options,
        layout=LinearMultiRecordOptions(),
        output=RenderOutputRequest(output_prefix="file", output_directory=tmp_path, overwrite=True),
    )
    render_request(request)
    planned = _materialized_record_inputs(resolve_request_records(request))
    path.unlink()
    live, replay = _replayed(
        LinearDiagramRequest(records=planned, options=options, layout=LinearMultiRecordOptions(),
                             output=_output(tmp_path)),
        tmp_path,
    )
    assert live == (tmp_path / "file.svg").read_text()
    assert replay == live


def _two_windows(source, tmp_path: Path, second: str) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=(
            RecordInput(
                source=source,
                windows=(RecordWindow(5001, 15000, False), RecordWindow.from_region(parse_region_specs([second])[0])),
                record_key="w",
            ),
        ),
        options=LinearDiagramOptions(config_overrides=_TRACKS),
        layout=LinearMultiRecordOptions(),
        output=_output(tmp_path),
    )


def _saved_windows(tmp_path: Path) -> list[tuple[str, int]]:
    document = read_session_container(tmp_path / "saved.gbdraw").document
    return [(item["recordKey"], len(item.get("windows", ()))) for item in document["drawings"][0]["renderRequest"]["records"]]


@pytest.mark.parametrize("second", ["src:40001-52000", "src:40001-52000:rc"], ids=["forward", "reversed"])
def test_a_session_keeps_the_windows_of_an_in_memory_record_as_one_record(second: str, tmp_path: Path) -> None:
    """Two windows of one in-memory record, given as the planner's units, are written as that record with both
    windows, so the replay draws the gap between them and each window's flanks (OV-2923, D-64)."""

    request = resolve_request(_two_windows(InMemoryRecordSource(_record("src", "linear")), tmp_path, second))
    live, replay = _replayed(request, tmp_path)
    assert live.count('data-gbdraw-row-mark="gap"') == 1
    assert _saved_windows(tmp_path) == [("w", 2)]
    assert replay == live


@pytest.mark.parametrize("second", ["src:40001-52000", "src:40001-52000:rc"], ids=["forward", "reversed"])
def test_a_session_keeps_the_windows_of_a_record_whose_file_is_gone(second: str, tmp_path: Path) -> None:
    """The planner's units of a file record with two windows, saved after the file is gone, are written as the
    record they were cut from with both windows (OV-2923 remaining work)."""

    path = tmp_path / "one.gb"
    SeqIO.write([_record("src", "linear")], path, "genbank")
    resolved = resolve_request(_two_windows(GenBankInputSource(str(path)), tmp_path, second))
    path.unlink()
    live, replay = _replayed(resolved, tmp_path)
    assert live.count('data-gbdraw-row-mark="gap"') == 1
    assert _saved_windows(tmp_path) == [("w", 2)]
    assert replay == live


_MITO = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"


def _feature_ids(svg: str) -> int:
    return len(set(re.findall(r'data-gbdraw-feature-id="([^"]+)"', svg)))


def _edit(record: SeqRecord) -> None:
    """A user's edit: one new CDS comes in beside the features the record has, and the record is renamed."""

    record.features.append(
        SeqFeature(FeatureLocation(10, 400, 1), type="CDS", qualifiers={"product": ["USEREDIT"], "gene": ["userEdit"]})
    )
    record.id = "edited"


def _resolved_mito(tmp_path: Path, windows: tuple[RecordWindow, ...] = ()) -> SeqRecord:
    resolved = resolve_request(LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(str(_MITO)), record_key="mito", windows=windows),),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
    ))
    source = resolved.records[0].source
    assert isinstance(source, InMemoryRecordSource) and held_record_origin(source.record) is not None
    return source.record


@pytest.mark.parametrize("windows", [(), (RecordWindow(100, 900),)], ids=["whole", "window"])
def test_a_session_saves_an_edited_copy_of_a_planned_record_as_itself(
    windows: tuple[RecordWindow, ...], tmp_path: Path
) -> None:
    """A planned record, copied, edited and given again, is saved as the edited record, not as its file (OV-2923)."""

    import copy

    record = copy.deepcopy(_resolved_mito(tmp_path, windows))
    _edit(record)
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(record), record_key="edited"),),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    assert _feature_ids(live) == _feature_ids(replay) > 1  # the new CDS beside the record's features
    # An edited crop is saved as the record it was cut from with its edits, and its window (OV-3034).
    assert replay == live


def test_an_edited_copy_of_a_planned_record_names_no_file(tmp_path: Path) -> None:
    """A planned record names its file when given again unchanged, and none once edited (losat_source_ids)."""

    import copy

    from gbdraw.comparisons.losat_jobs import record_source_paths

    def replanned(record: SeqRecord) -> SeqRecord:
        resolved = resolve_request(LinearDiagramRequest(
            records=(RecordInput(InMemoryRecordSource(record), record_key="again"),),
            options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
        ))
        return resolved.records[0].source.record  # type: ignore[union-attr]

    planned = _resolved_mito(tmp_path)
    assert record_source_paths(replanned(copy.deepcopy(planned))) == (str(_MITO),)
    edited = copy.deepcopy(planned)
    _edit(edited)
    assert record_source_paths(replanned(edited)) == ()


def test_a_session_saves_a_planned_record_edited_in_place_as_itself(tmp_path: Path) -> None:
    """A planned record edited in place, its projected source kept, no longer reads as its file's record."""

    resolved = resolve_request(LinearDiagramRequest(
        records=(RecordInput(GenBankInputSource(str(_MITO)), record_key="mito"),),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
    ))
    source = resolved.records[0].source
    assert isinstance(source, InMemoryRecordSource)
    _edit(source.record)
    live, replay = _replayed(resolved, tmp_path)
    assert _feature_ids(live) > 1
    assert replay == live


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
def test_a_session_redraws_an_edited_crop_from_its_source_with_the_edits(spec: str, tmp_path: Path) -> None:
    """A crop whose features and name were edited after its cut is saved as the record it was cut from, with the
    edits mapped into the crop's span (reversed back for ``:rc``), and the crop's window, so the replay draws the
    edited crop with its flanks and its source coordinates (OV-3034)."""

    crop = apply_region_specs([_record("src", "linear")], parse_region_specs([spec]))[0]
    _edit(crop)
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(crop), record_key="r"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    assert _feature_ids(live) == 2  # the record's CDS and the new one
    [saved] = _saved_records(tmp_path)
    assert saved["windows"] == [
        {"start": 20001, "end": 40000, "reverseComplement": spec.endswith(":rc"), "column": None}
    ]
    assert replay == live


def test_an_edited_reverse_complemented_crop_of_an_rna_record_is_saved_on_its_source(tmp_path: Path) -> None:
    """A reverse-complemented crop of an RNA record is complemented as RNA; its edits are written on its source,
    whose window is complemented as the crop was (OV-3332: it was compared as DNA and saved as its own bytes)."""

    record = _record("src", "linear")
    record.seq = Seq(str(record.seq).replace("T", "U"))
    record.annotations["molecule_type"] = "ss-RNA"
    crop = apply_region_specs([record], parse_region_specs(["src:20001-40000:rc"]))[0]
    assert "U" in str(crop.seq)
    _edit(crop)
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(crop), record_key="r"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    [saved] = _saved_records(tmp_path)
    assert saved["windows"] == [{"start": 20001, "end": 40000, "reverseComplement": True, "column": None}]
    assert replay == live


def _flanked_source() -> SeqRecord:
    """A circular 60 kb record with a CDS across base 20,000, a CDS joined across the origin, and one inside 20-40 kb."""

    record = _record("src", "circular")
    record.features = [
        SeqFeature(FeatureLocation(0, 60_000, 1), type="source", qualifiers={"organism": ["test"]}),
        SeqFeature(FeatureLocation(19_000, 21_000, 1), type="CDS", qualifiers={"locus_tag": ["S1"], "product": ["s"]}),
        SeqFeature(FeatureLocation(59_000, 60_000, 1) + FeatureLocation(0, 500, 1), type="CDS",
                   qualifiers={"locus_tag": ["ORI"], "product": ["o"]}),
        SeqFeature(FeatureLocation(25_000, 26_000, 1), type="CDS", qualifiers={"locus_tag": ["I1"], "product": ["i"]}),
    ]
    return record


def _saved_source(tmp_path: Path) -> SeqRecord:
    import io

    container = read_session_container(tmp_path / "saved.gbdraw")
    [saved] = container.document["drawings"][0]["renderRequest"]["records"]
    with container.resource(saved["source"]["resourceId"]).open_bytes() as stream:
        data = stream.read()
    return SeqIO.read(io.StringIO(data.decode("utf-8")), "genbank")


def _locus(feature: SeqFeature) -> str:
    return feature.qualifiers.get("locus_tag", [feature.type])[0]


_CROP_EDITS = {
    # Only a qualifier of a CDS inside the crop changes.
    "product": lambda crop: next(f for f in crop.features if _locus(f) == "I1").qualifiers.update(product=["edited"]),
    # Only the crop's name changes.
    "rename": lambda crop: setattr(crop, "id", "renamed"),
    # A CDS inside the crop is taken out.
    "remove": lambda crop: setattr(crop, "features", [f for f in crop.features if _locus(f) != "I1"]),
}


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("edit", sorted(_CROP_EDITS))
def test_an_edited_crop_is_saved_with_only_its_edits_on_its_source(
    mode: str, spec: str, edit: str, tmp_path: Path
) -> None:
    """The source of an edited crop keeps the features the crop cut at its ends whole, and the features outside the
    crop (one joined across the origin too); only what the edit changed is written to it, so the replay draws the
    live figure (OV-3034)."""

    import copy

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    _CROP_EDITS[edit](crop)
    record = RecordInput(source=InMemoryRecordSource(crop), record_key="r")
    request = (
        LinearDiagramRequest(records=(record,), options=LinearDiagramOptions(config_overrides=_TRACKS),
                             output=_output(tmp_path))
        if mode == "linear"
        else CircularDiagramRequest(records=(record,), options=CircularDiagramOptions(config_overrides=_TRACKS),
                                    output=_output(tmp_path))
    )
    live, replay = _replayed(request, tmp_path)
    assert replay == live
    [saved] = _saved_records(tmp_path)
    assert [(window["start"], window["end"]) for window in saved["windows"]] == [(20001, 40000)]
    source = _saved_source(tmp_path)
    original = copy.deepcopy(_flanked_source())
    expected = [(_locus(feature), str(feature.location), feature.qualifiers.get("product"))
                for feature in original.features]
    if edit == "product":
        expected = [(locus, location, ["edited"] if locus == "I1" else product) for locus, location, product in expected]
    if edit == "remove":
        expected = [item for item in expected if item[0] != "I1"]
    # In the source's order: the order of the saved source gives each feature its source feature index on replay.
    assert [(_locus(feature), str(feature.location), feature.qualifiers.get("product")) for feature in source.features] \
        == expected
    assert source.id == ("renamed" if edit == "rename" else "src")


def _request(mode: str, crop: SeqRecord, tmp_path: Path) -> LinearDiagramRequest | CircularDiagramRequest:
    record = RecordInput(source=InMemoryRecordSource(crop), record_key="r")
    if mode == "linear":
        return LinearDiagramRequest(records=(record,), options=LinearDiagramOptions(config_overrides=_TRACKS),
                                    output=_output(tmp_path))
    return CircularDiagramRequest(records=(record,), options=CircularDiagramOptions(config_overrides=_TRACKS),
                                  output=_output(tmp_path))


def _drawn_ids(svg: str) -> set[str]:
    return set(re.findall(r'data-gbdraw-feature-id="([^"]+)"', svg))


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("place", [3, 1], ids=["appended", "inserted"])
def test_a_feature_added_to_a_crop_that_keeps_others_draws_and_replays(
    mode: str, spec: str, place: int, tmp_path: Path
) -> None:
    """A feature added to a crop has no source feature: it takes a source feature index no feature cut from the
    source can hold, so the live figure draws it beside the cut features and the replay draws the live figure, in
    the same order (OV-3141, OV-3229). Appended, its list position (3) is the index of the CDS the crop cut from the
    source's fourth feature."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    crop.features.insert(
        place,
        SeqFeature(FeatureLocation(1000, 2000, 1), type="CDS", qualifiers={"locus_tag": ["NEW"], "product": ["n"]}),
    )
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert replay == live
    assert len(_drawn_ids(live)) == 3  # S1, I1 and the added CDS (the source feature is not drawn)
    assert sorted(_locus(feature) for feature in _saved_source(tmp_path).features) == \
        ["I1", "NEW", "ORI", "S1", "source"]


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
def test_a_feature_moved_inside_a_crop_has_one_id_live_and_on_replay(mode: str, spec: str, tmp_path: Path) -> None:
    """A feature moved inside a crop is a feature at its new place, as a moved feature of a whole record is: the live
    figure and the replay give it the ID of its new source location (OV-3142)."""

    from gbdraw.features.ids import compute_feature_hash_from_location_parts

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    first = tmp_path / "unedited"
    first.mkdir()
    unedited = _drawn_ids(_replayed(_request(mode, crop, first), first)[0])
    moved = next(feature for feature in crop.features if _locus(feature) == "I1")
    moved.location = FeatureLocation(7000, 8000, moved.location.strand)
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert replay == live
    # The crop's 7,001-8,000 is the source's 27,001-28,000, or 32,001-33,000 in a reversed crop (the CDS stays on the
    # source's forward strand).
    start = 32_000 if spec.endswith(":rc") else 27_000
    expected = compute_feature_hash_from_location_parts("CDS", [(start, start + 1000, 1)], record_id="src")
    drawn = _drawn_ids(live)
    assert expected in drawn
    assert len(drawn) == len(unedited) and len(drawn - unedited) == 1


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
def test_a_moved_feature_reads_its_source_location_from_its_new_place(spec: str) -> None:
    """Once a crop's feature is moved, its source location parts and source strand are those of its new place; a
    record gbdraw loads them from keeps them for every later reader (OV-3142)."""

    import copy

    from gbdraw.core.record_metadata import _source_feature_anchor_profile, _source_feature_location_parts
    from gbdraw.features.source import capture_source_features

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    moved = next(feature for feature in crop.features if _locus(feature) == "I1")
    moved.location = FeatureLocation(7000, 8000, moved.location.strand)
    planned = copy.deepcopy(crop)
    capture_source_features(planned)
    [feature] = [feature for feature in planned.features if _locus(feature) == "I1"]
    start = 32_000 if spec.endswith(":rc") else 27_000
    assert _source_feature_location_parts(feature) == ((start, start + 1000, 1),)
    assert _source_feature_anchor_profile(feature).strand == "+"


def test_a_copy_of_a_cut_feature_takes_a_source_feature_index_of_its_own(tmp_path: Path) -> None:
    """A copy of a feature the crop cut carries that feature's source feature index; the copy takes the next free
    index, so the live figure draws both (OV-3141)."""

    import copy

    crop = apply_region_specs([_flanked_source()], parse_region_specs(["src:20001-40000"]))[0]
    crop.features.append(copy.deepcopy(next(feature for feature in crop.features if _locus(feature) == "I1")))
    render_request(_request("linear", crop, tmp_path))
    # The two copies of I1 are drawn as two instances of one feature, named by their source feature indexes.
    assert sorted(set(re.findall(r'id="f[0-9a-f]+__instance_(\d+)_', (tmp_path / "live.svg").read_text()))) == ["3", "4"]


def _record_facts(svg: str) -> tuple[str, list[str], list[str], set[str]]:
    """What a replay that saved the crop as its own bytes changed: the figure width, the record's subtitle and source
    span attributes, and the feature IDs."""

    width = re.search(r'<svg[^>]*\swidth="([^"]+)"', svg)
    assert width is not None
    return (
        width.group(1),
        re.findall(r'>([0-9,]+(?: bp|-[0-9,]+))</text>', svg),
        re.findall(r'data-gbdraw-record-source-(?:start|end|step)="[^"]*"', svg),
        _drawn_ids(svg),
    )


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("copied", ["I1", "S1"], ids=["inside", "cut-at-start"])
def test_a_copy_of_a_cut_feature_is_saved_as_an_added_feature_and_replays_as_drawn(
    mode: str, spec: str, copied: str, tmp_path: Path
) -> None:
    """A copy of a feature the crop cut is a feature added at the copy's place, live and in the saved Session: the
    crop is saved as its source with the copy added, so the replay keeps the crop's subtitle, width, source span
    and feature IDs (OV-3247). A copy of a feature the crop cut at its start is the part inside the crop."""

    import copy

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    crop.features.append(copy.deepcopy(next(feature for feature in crop.features if _locus(feature) == copied)))
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert _record_facts(replay) == _record_facts(live)
    assert replay == live
    [saved] = _saved_records(tmp_path)
    assert [(window["start"], window["end"]) for window in saved["windows"]] == [(20001, 40000)]
    assert sorted(_locus(feature) for feature in _saved_source(tmp_path).features) == \
        sorted(["source", "S1", "ORI", "I1", copied])



@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
def test_a_feature_copied_from_another_crop_of_the_source_is_an_added_feature(
    mode: str, spec: str, tmp_path: Path
) -> None:
    """A feature copied from another crop of the same record carries a source feature index this crop did not cut;
    it is a feature added at its place, live and in the saved Session, so the replay keeps the crop's subtitle,
    width, source span, feature IDs and paint order (OV-3256)."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    other = apply_region_specs([_flanked_source()], parse_region_specs(["src:1-10000"]))[0]
    crop.features.append(copy.deepcopy(next(feature for feature in other.features if _locus(feature) == "ORI")))
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert _record_facts(replay) == _record_facts(live)
    assert replay == live
    [saved] = _saved_records(tmp_path)
    assert [(window["start"], window["end"]) for window in saved["windows"]] == [(20001, 40000)]
    # The source keeps ORI outside the window and gains the copy as an added feature.
    assert sorted(_locus(feature) for feature in _saved_source(tmp_path).features) == \
        sorted(["source", "S1", "ORI", "I1", "ORI"])


_ORIGIN_EDITS = {
    # Only a qualifier of the CDS joined across the origin changes.
    "product": lambda crop: next(f for f in crop.features if _locus(f) == "ORI").qualifiers.update(product=["edited"]),
    # A copy of the CDS joined across the origin.
    "copy": lambda crop: crop.features.append(
        copy.deepcopy(next(f for f in crop.features if _locus(f) == "ORI"))
    ),
    # A CDS added across the source's origin (crop bases 9,501-10,500 are source bases 59,501-500).
    "added": lambda crop: crop.features.append(
        SeqFeature(FeatureLocation(9_500, 10_500, 1), type="CDS", qualifiers={"locus_tag": ["NEW"], "product": ["n"]})
    ),
    # A partial CDS added across the source's origin, open at both ends: its partial ends stay on replay.
    "partial": lambda crop: crop.features.append(
        SeqFeature(FeatureLocation(BeforePosition(9_900), AfterPosition(10_200), 1), type="CDS",
                   qualifiers={"locus_tag": ["NEW"], "product": ["n"]})
    ),
}


@pytest.mark.parametrize(
    ("location", "expected"),
    [
        (FeatureLocation(BeforePosition(59_900), AfterPosition(60_200), 1), "join{[<59900:60000](+), [0:>200](+)}"),
        (FeatureLocation(BeforePosition(59_900), AfterPosition(60_200), -1), "join{[0:>200](-), [<59900:60000](-)}"),
        (FeatureLocation(BeforePosition(60_100), AfterPosition(60_400), 1), "[<100:>400](+)"),
        (FeatureLocation(BeforePosition(100), AfterPosition(900), -1), "[<100:>900](-)"),
    ],
    ids=["across-plus", "across-minus", "past", "before"],
)
def test_a_place_read_past_the_origin_is_wrapped_onto_the_circle_with_its_partial_ends(location, expected) -> None:
    """A place read past the origin of a circle goes back onto it; a part across the origin is split there, the cut
    points exact and the outer ends keeping their partial positions, and reading it round the origin gives the place
    back (OV-3262)."""

    from gbdraw.io.regions import _unwrapped_location, _wrapped_location

    wrapped = _wrapped_location(location, 60_000)
    assert str(wrapped) == expected
    if isinstance(wrapped, CompoundLocation):
        assert str(_unwrapped_location(wrapped, 60_000)) == str(location)


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize(
    "spec", ["src:50001-10000:origin", "src:50001-10000:origin:rc"], ids=["forward", "reversed"]
)
@pytest.mark.parametrize("edit", sorted(_ORIGIN_EDITS))
def test_an_edited_crop_across_the_origin_is_saved_with_its_edits_on_its_source(
    mode: str, spec: str, edit: str, tmp_path: Path
) -> None:
    """An edited crop across the origin is saved as its source with the edits and its window, as a crop that does
    not cross the origin is, so the replay keeps its subtitle, width, source span and feature IDs (OV-3262)."""

    from gbdraw.io.regions import crop_record_to_region

    crop = crop_record_to_region(_flanked_source(), parse_region_specs([spec])[0], circular=True)
    _ORIGIN_EDITS[edit](crop)
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert _record_facts(replay) == _record_facts(live)
    assert replay == live
    [saved] = _saved_records(tmp_path)
    assert [(window["start"], window["end"]) for window in saved["windows"]] == [(50001, 10000)]


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize(
    "inner", ["src:5001-15000:rc", "src:5001-15000", "src:12001-15000", "src:12001-15000:rc"],
    ids=["across-reversed", "across-forward", "past-forward", "past-reversed"],
)
def test_a_crop_of_a_crop_across_the_origin_draws_and_replays(mode: str, inner: str, tmp_path: Path) -> None:
    """A crop of a crop across the origin reads its source bases round the origin, whether it runs across the origin
    or lies past it, forward or reverse-complemented: it draws, and its Session replays it as drawn (OV-3270)."""

    from gbdraw.io.regions import crop_record_to_region

    outer = crop_record_to_region(_flanked_source(), parse_region_specs(["src:50001-10000:origin"])[0], circular=True)
    crop = crop_record_to_region(outer, parse_region_specs([inner])[0])
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert replay == live


def _out_of_start_order_source() -> SeqRecord:
    """A circular 60 kb record whose overlapping features are listed out of start order, inside 20-40 kb and across
    the origin."""

    record = _record("src", "circular")
    record.features = [
        SeqFeature(FeatureLocation(0, 60_000, 1), type="source", qualifiers={"organism": ["test"]}),
        SeqFeature(FeatureLocation(25_000, 27_000, 1), type="CDS", qualifiers={"locus_tag": ["LATE"]}),
        SeqFeature(FeatureLocation(24_000, 26_000, 1), type="tRNA", qualifiers={"locus_tag": ["EARLY"]}),
        SeqFeature(FeatureLocation(59_000, 60_000, 1) + FeatureLocation(0, 500, 1), type="CDS",
                   qualifiers={"locus_tag": ["ORI"]}),
        SeqFeature(FeatureLocation(58_000, 60_000, 1) + FeatureLocation(0, 1_000, 1), type="tRNA",
                   qualifiers={"locus_tag": ["ORI2"]}),
    ]
    return record


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize(
    ("outer", "inner"),
    [("src:10001-50000:rc", "src:1001-30000:rc"), ("src:50001-10000:origin:rc", "src:5001-15000:rc")],
    ids=["inside", "across-the-origin"],
)
def test_a_reverse_complemented_crop_of_a_reverse_complemented_crop_draws_in_source_order(
    mode: str, outer: str, inner: str, tmp_path: Path
) -> None:
    """Reversed twice, a crop is back in its source's orientation and lists its features in source order, as the
    forward cut its Session replays does, so the replay paints overlapping features as the live figure does
    (OV-3271)."""

    from gbdraw.io.regions import crop_record_to_region

    crop = crop_record_to_region(
        crop_record_to_region(_out_of_start_order_source(), parse_region_specs([outer])[0], circular=True),
        parse_region_specs([inner])[0],
    )
    expected = ["source", "LATE", "EARLY"] if "10001" in outer else ["source", "ORI", "ORI2"]
    assert [_locus(feature) for feature in crop.features] == expected
    live, replay = _replayed(_request(mode, crop, tmp_path), tmp_path)
    assert replay == live

def _twins_source() -> SeqRecord:
    """A linear 60 kb record with a CDS X and two identical CDS B and C (same type and location) inside 20-40 kb."""

    record = _record("src", "linear")
    record.features = [
        SeqFeature(FeatureLocation(0, 60_000, 1), type="source", qualifiers={"organism": ["test"]}),
        SeqFeature(FeatureLocation(21_000, 22_000, 1), type="CDS", qualifiers={"locus_tag": ["X"], "product": ["x"]}),
        SeqFeature(FeatureLocation(30_000, 31_000, 1), type="CDS", qualifiers={"locus_tag": ["B"], "product": ["b"]}),
        SeqFeature(FeatureLocation(30_000, 31_000, 1), type="CDS", qualifiers={"locus_tag": ["C"], "product": ["c"]}),
    ]
    return record


def _instances(svg: str) -> list[str]:
    """The rendered IDs of identical features, named by their source feature indexes (``__instance_<index>``)."""
    return sorted(set(re.findall(r'\bid="(f[0-9a-f]+__instance_\d+)_', svg)))


_TWIN_EDITS = {
    "rename": lambda crop: setattr(crop, "id", "renamed"),
    # A CDS inserted before the crop's features.
    "insert": lambda crop: crop.features.insert(
        1, SeqFeature(FeatureLocation(1000, 2000, 1), type="CDS", qualifiers={"locus_tag": ["NEW"], "product": ["n"]})
    ),
    # Three CDS appended out of start order, two of them identical.
    "append": lambda crop: crop.features.extend(
        SeqFeature(FeatureLocation(start, start + 1000, 1), type="CDS", qualifiers={"locus_tag": [tag]})
        for start, tag in ((5000, "P"), (1000, "Q1"), (1000, "Q2"))
    ),
}


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("edit", sorted(_TWIN_EDITS))
def test_identical_features_of_an_edited_crop_keep_their_ids_on_replay(spec: str, edit: str, tmp_path: Path) -> None:
    """Identical features are told apart by their source feature indexes; an edit that removes nothing keeps every
    index on replay, so the replay names them as the live figure does: the saved source keeps every kept feature at
    its source place and lists the added ones after them."""

    crop = apply_region_specs([_twins_source()], parse_region_specs([spec]))[0]
    _TWIN_EDITS[edit](crop)
    live, replay = _replayed(_request("linear", crop, tmp_path), tmp_path)
    # B and C at least; in a forward crop a CDS added at 1,001-2,000 lies where X lies, so X and it are identical too.
    assert len(_instances(live)) >= 2
    assert _instances(replay) == _instances(live)
    assert replay == live


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
def test_a_feature_hidden_in_an_edited_crop_stays_hidden_after_save_and_load(spec: str, tmp_path: Path) -> None:
    """A feature override names one of two identical features by its source feature index; the replay of the
    saved Session hides the feature the live figure hid."""

    from gbdraw.features.overrides import FeatureOverride

    crop = apply_region_specs([_twins_source()], parse_region_specs([spec]))[0]
    crop.id = "renamed"
    probe = tmp_path / "probe"
    probe.mkdir()
    render_request(_request("linear", crop, probe))
    [feature_hash] = {instance.split("__")[0] for instance in _instances((probe / "live.svg").read_text())}
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(crop), record_key="r"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS, feature_overrides=(FeatureOverride(
            record_key="r", biological_feature_id=f"{feature_hash}~3", feature_visibility="off"),)),
        output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    assert _instances(live) == [f"{feature_hash}__instance_2"]  # C, source feature 3, is hidden
    assert replay == live


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("edited", [False, True], ids=["alone", "with-an-edit"])
def test_a_cut_feature_given_an_equal_location_keeps_its_id(spec: str, edited: bool, tmp_path: Path) -> None:
    """A feature whose location is replaced by an equal one has not moved, live as in a Session: the feature the crop
    cut at its start keeps the ID of its whole source location in both, alone or beside an edit that makes the
    Session write the crop's source with its edits: one rule, a location equal as GenBank writes it, says whether a
    feature moved."""

    import copy

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    first = tmp_path / "unedited"
    first.mkdir()
    unedited = _drawn_ids(_replayed(_request("linear", crop, first), first)[0])
    cut = next(feature for feature in crop.features if _locus(feature) == "S1")
    cut.location = copy.deepcopy(cut.location)
    if edited:
        next(feature for feature in crop.features if _locus(feature) == "I1").qualifiers["note"] = ["edited"]
    live, replay = _replayed(_request("linear", crop, tmp_path), tmp_path)
    assert _drawn_ids(live) == unedited
    assert replay == live


def test_the_web_view_parts_of_a_crop_read_the_source_feature_index_owner() -> None:
    """The Web helper that indexes view location parts by source feature index reads the index owner, so a crop with
    an added feature gives each feature its own parts instead of two features one index."""

    from gbdraw.core.record_metadata import _feature_source_index_map

    helpers_js = (Path(__file__).parents[1] / "gbdraw" / "web" / "js" / "app" / "python-helpers.js").read_text(
        encoding="utf-8"
    )
    namespace: dict[str, object] = {}
    exec(helpers_js.split("`", 1)[1].rsplit("`", 1)[0], namespace)
    crop = apply_region_specs([_flanked_source()], parse_region_specs(["src:20001-40000"]))[0]
    crop.features.insert(1, SeqFeature(FeatureLocation(1000, 2000, 1), type="CDS", qualifiers={"locus_tag": ["NEW"]}))
    indexes = _feature_source_index_map(crop.features)
    expected = {
        indexes[id(feature)]: tuple((int(part.start), int(part.end), part.strand) for part in feature.location.parts)
        for feature in crop.features
    }
    assert namespace["_web_feature_view_hash_parts_index"](crop) == expected  # type: ignore[operator]


def test_an_edited_reverse_complemented_crop_draws_a_feature_with_no_location(tmp_path: Path) -> None:
    """A feature with no location is left out of the drawing, as in any record; ordering an edited crop's features for
    its replay sorts it last instead of stopping the render."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs(["src:20001-40000:rc"]))[0]
    crop.features.append(SeqFeature(None, type="misc_feature", qualifiers={"note": ["no place"]}))
    render_request(_request("linear", crop, tmp_path))
    assert len(_drawn_ids((tmp_path / "live.svg").read_text())) == 2


def test_an_edited_crop_saved_as_its_own_bytes_lists_its_features_as_drawn(tmp_path: Path) -> None:
    """A crop saved as its own bytes lists its features in the order the live figure draws them, so the replay paints
    them in that order too."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs(["src:20001-40000"]))[0]
    # A crop whose bases were edited is saved as its own bytes.
    crop.seq = Seq("A" * 100) + crop.seq[100:]
    crop.features.insert(
        1, SeqFeature(FeatureLocation(1000, 2000, 1), type="CDS", qualifiers={"locus_tag": ["NEW"], "product": ["n"]})
    )
    request = _request("linear", crop, tmp_path)
    planned = resolve_request(request).records[0].source.record  # type: ignore[union-attr]
    _replayed(request, tmp_path)
    [saved] = _saved_records(tmp_path)
    assert not saved.get("windows")
    assert [_locus(feature) for feature in _saved_source(tmp_path).features] == \
        [_locus(feature) for feature in planned.features]


_END_CUT_EDITS = {
    # S1, which the crop cut at its start, is moved inside the crop.
    "moved": lambda feature: setattr(feature, "location", FeatureLocation(100, 900, feature.location.strand)),
    # S1 changes its type only.
    "retyped": lambda feature: setattr(feature, "type", "gene"),
}


@pytest.mark.parametrize("spec", ["src:20001-40000", "src:20001-40000:rc"], ids=["forward", "reversed"])
@pytest.mark.parametrize("edit", sorted(_END_CUT_EDITS))
def test_a_session_redraws_an_edited_feature_the_crop_cut_at_its_end(spec: str, edit: str, tmp_path: Path) -> None:
    """A feature the crop cut at an end, once moved, is a feature at its new place; once retyped, it keeps its whole
    source location. Either way the crop is saved as its source with the edit, and the replay draws the live figure
    with the same feature IDs (OV-3246)."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs([spec]))[0]
    _END_CUT_EDITS[edit](next(feature for feature in crop.features if _locus(feature) == "S1"))
    live, replay = _replayed(_request("linear", crop, tmp_path), tmp_path)
    assert _drawn_ids(replay) == _drawn_ids(live)
    assert replay == live
    [saved] = _saved_records(tmp_path)
    assert [(window["start"], window["end"]) for window in saved["windows"]] == [(20001, 40000)]


def test_a_session_saves_a_crop_whose_bases_were_edited_as_itself(tmp_path: Path) -> None:
    """A crop whose bases were edited is no longer a span of its source; it is saved as its own bytes."""

    crop = apply_region_specs([_record("src", "linear")], parse_region_specs(["src:20001-40000"]))[0]
    crop.seq = Seq("A" * 100) + crop.seq[100:]
    crop.features = []
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(crop), record_key="r"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), output=_output(tmp_path),
    )
    _live, _replay = _replayed(request, tmp_path)
    [saved] = _saved_records(tmp_path)
    assert not saved.get("windows")


@pytest.mark.parametrize("case", ["crop", "rc", "wcrop", "wrc"])
def test_a_session_redraws_a_planned_record_given_again_with_a_new_window_or_orientation(
    case: str, tmp_path: Path
) -> None:
    """A planned record given again with a window, or reverse-complemented, replays as drawn."""

    record = _resolved_mito(tmp_path, (RecordWindow(1000, 9000),) if case.startswith("w") else ())
    request = LinearDiagramRequest(
        records=(
            RecordInput(InMemoryRecordSource(record), record_key="m2", windows=(RecordWindow(2000, 3000),))
            if case.endswith("crop")
            else RecordInput(
                InMemoryRecordSource(record), record_key="m2", presentation=RecordPresentation(reverse_complement=True)
            ),
        ),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    assert replay == live


def test_a_session_redraws_a_crop_of_a_reverse_complemented_record(tmp_path: Path) -> None:
    """A crop of a planned, reverse-complemented record is saved as the record's forward form and the crop's window."""

    resolved = resolve_request(LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(_record("src", "linear")), record_key="w",
                             presentation=RecordPresentation(reverse_complement=True)),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), layout=LinearMultiRecordOptions(),
        output=_output(tmp_path),
    ))
    crop = apply_region_specs([resolved.records[0].source.record], parse_region_specs(["src:5001-15000"]))[0]  # type: ignore[union-attr]
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(crop), record_key="c"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), layout=LinearMultiRecordOptions(),
        output=_output(tmp_path),
    )
    live, replay = _replayed(request, tmp_path)
    assert _saved_windows(tmp_path) == [("c", 1)]
    assert replay == live


def test_an_edited_window_unit_is_no_longer_a_window(tmp_path: Path) -> None:
    """Of two planned windows of one record given again, the edited one is drawn and saved as a record of its own,
    so live and replay agree that no gap mark joins them."""

    import copy

    resolved = resolve_request(_two_windows(InMemoryRecordSource(_record("src", "linear")), tmp_path, "src:40001-52000"))
    first, second = (copy.deepcopy(record.source.record) for record in resolved.records)  # type: ignore[union-attr]
    second.features.append(SeqFeature(FeatureLocation(100, 400, 1), type="CDS", qualifiers={"product": ["EDIT"]}))
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(first), record_key="a"),
                 RecordInput(InMemoryRecordSource(second), record_key="b")),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
    )
    replanned = resolve_request(request)
    assert [record.source.record.annotations.get("gbdraw_window_of") for record in replanned.records] == ["w", None]  # type: ignore[union-attr]
    live, replay = _replayed(request, tmp_path)
    assert live.count('data-gbdraw-row-mark="gap"') == replay.count('data-gbdraw-row-mark="gap"') == 0


def test_a_file_qualified_region_selects_an_in_memory_record_read_from_that_file(tmp_path: Path) -> None:
    """A record read with ``read_genbank`` keeps the file stamps the loader wrote, so a region named by its file
    still selects it when the record is given in memory (as on dev)."""

    import gbdraw

    source = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    record = gbdraw.read_genbank(str(source))[0]
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(), layout=LinearMultiRecordOptions(), output=_output(tmp_path),
        record_options=RecordCollectionOptions(regions=tuple(parse_region_specs(["HmmtDNA.gbk:NC_012920.1:1000-5000"]))),
    )
    resolved = resolve_request(request)
    drawn = resolved.records[0].source.record  # type: ignore[union-attr]
    assert len(drawn.seq) == 4001
    assert drawn.annotations.get("gbdraw_source_basename") == "HmmtDNA.gbk"
