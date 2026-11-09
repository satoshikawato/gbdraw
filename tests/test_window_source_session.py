"""A crop's windows read the bases beside it from the record it was cut from (OV-2890); a Session keeps them.

A crop given in memory, or a crop of a file that is no longer at hand, is written to a Session as the record it was
cut from and the crop's window in it, so the replay crops it again and draws what the live figure drew (D-64).
"""

from __future__ import annotations

import json
import random
import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
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
    path = tmp_path / "saved.gbdraw-session.json"
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
    document = json.loads((tmp_path / "saved.gbdraw-session.json").read_text())
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
    path = tmp_path / "saved.gbdraw-session.json"
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
    document = json.loads((tmp_path / "saved.gbdraw-session.json").read_text())
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
    """A user's edit: every feature but the source goes, one new CDS comes in, and the record is renamed."""

    record.features = [feature for feature in record.features if feature.type == "source"]
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
    assert _feature_ids(live) == _feature_ids(replay) == 1  # the new CDS alone
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
    assert _feature_ids(live) == 1
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
    assert _feature_ids(live) == 1  # the new CDS alone
    [saved] = _saved_records(tmp_path)
    assert saved["windows"] == [
        {"start": 20001, "end": 40000, "reverseComplement": spec.endswith(":rc"), "column": None}
    ]
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
    import base64
    import io

    document = json.loads((tmp_path / "saved.gbdraw-session.json").read_text())
    [saved] = document["drawings"][0]["renderRequest"]["records"]
    data = base64.b64decode(document["resources"][saved["source"]["resourceId"]]["data"])
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
    expected = {_locus(feature): (str(feature.location), feature.qualifiers.get("product")) for feature in original.features}
    if edit == "product":
        expected["I1"] = (expected["I1"][0], ["edited"])
    if edit == "remove":
        del expected["I1"]
    assert {_locus(feature): (str(feature.location), feature.qualifiers.get("product")) for feature in source.features} \
        == expected
    assert source.id == ("renamed" if edit == "rename" else "src")


def test_a_session_saves_a_crop_that_moved_a_feature_cut_at_its_end_as_itself(tmp_path: Path) -> None:
    """A feature the crop cut at an end has no place on the source once moved; the crop is saved as its own bytes."""

    crop = apply_region_specs([_flanked_source()], parse_region_specs(["src:20001-40000"]))[0]
    cut = next(feature for feature in crop.features if _locus(feature) == "S1")
    cut.location = FeatureLocation(100, 900, 1)
    request = LinearDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(crop), record_key="r"),),
        options=LinearDiagramOptions(config_overrides=_TRACKS), output=_output(tmp_path),
    )
    _replayed(request, tmp_path)
    [saved] = _saved_records(tmp_path)
    assert not saved.get("windows")
    assert [str(feature.location) for feature in _saved_source(tmp_path).features if _locus(feature) == "S1"] \
        == ["[100:900](+)"]


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
