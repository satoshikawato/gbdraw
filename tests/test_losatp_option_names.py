"""LOSATP option names: retired CLI flags, renamed Python fields, wire names.

Design: docs/internal/LOSAT_CLI_API_DESIGN_PROPOSAL_2026-10-03.md sections 3.2,
3.4, and 3.5 (Owner decisions D4, D5, D6, D8, D13, D16).
"""

from __future__ import annotations

import copy
import gzip
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from svgwrite import Drawing

import gbdraw.api.request_render as request_render_module
import gbdraw.linear as linear_cli_module
from gbdraw.api.requests import LinearDiagramRequest


FIXTURES = Path(__file__).parent / "fixtures" / "sessions"
GALLERY_SESSIONS = Path(__file__).parents[1] / "gbdraw" / "web" / "gallery" / "sessions"

# Design 3.5: every retired CLI input and the replacement its rejection names.
RETIRED_FLAGS = (
    (["--protein_blastp_mode", "none"], "omit"),
    (["--protein_blastp_mode", "pairwise"], "--losat losatp --losatp_mode pairwise"),
    (
        ["--protein_blastp_mode", "orthogroup"],
        "--losat losatp --losatp_mode similarity_groups",
    ),
    (["--protein_blastp_mode", "collinear"], "--losat losatp --losatp_mode collinear"),
    (["--losatp_bin", "custom-losat"], "--losat_bin"),
    (["--ncbi_blastp_bin", "/opt/blastp"], "--ncbi_blast_bin"),
    (["--losatp_threads", "2"], "--losat_threads"),
    (["--protein_blastp_max_hits", "3"], "--losatp_max_hits"),
    (["--protein_blastp_candidate_limit", "4"], "--losatp_max_target_seqs"),
    (["--align_orthogroup_feature", "og_1"], "--similarity_alignment_feature"),
    (["--protein_blastp_output", "raw.tsv"], "--losat_output_dir"),
)

# Design 3.5: legacy session argv rewrites (one token may become several).
LEGACY_REWRITES = (
    (["--protein_blastp_mode", "none"], []),
    (["--protein_blastp_mode", "pairwise"], ["--losat", "losatp", "--losatp_mode", "pairwise"]),
    (
        ["--protein_blastp_mode", "orthogroup"],
        ["--losat", "losatp", "--losatp_mode", "similarity_groups"],
    ),
    (
        ["--protein_blastp_mode=collinear"],
        ["--losat", "losatp", "--losatp_mode", "collinear"],
    ),
    (
        ["--protein-blastp-mode", "orthogroup"],
        ["--losat", "losatp", "--losatp_mode", "similarity_groups"],
    ),
    (["--losatp_bin", "x"], ["--losat_bin", "x"]),
    (["--losatp-bin", "x"], ["--losat_bin", "x"]),
    (["--losatp_bin=x"], ["--losat_bin=x"]),
    (["--ncbi_blastp_bin", "y"], ["--ncbi_blast_bin", "y"]),
    (["--ncbi-blastp-bin", "y"], ["--ncbi_blast_bin", "y"]),
    (["--losatp_threads", "2"], ["--losat_threads", "2"]),
    (["--losatp-threads", "2"], ["--losat_threads", "2"]),
    (["--protein_blastp_max_hits", "3"], ["--losatp_max_hits", "3"]),
    (["--protein-blastp-max-hits", "3"], ["--losatp_max_hits", "3"]),
    (["--protein_blastp_candidate_limit", "4"], ["--losatp_max_target_seqs", "4"]),
    (["--protein-blastp-candidate-limit", "4"], ["--losatp_max_target_seqs", "4"]),
    (["--align_orthogroup_feature", "og_1"], ["--similarity_alignment_feature", "og_1"]),
    (["--align-orthogroup-feature", "og_1"], ["--similarity_alignment_feature", "og_1"]),
)

# Old argv (as a 0.13.0 session stores it) and the new argv it must equal.
EQUIVALENT_ARGV = (
    (
        ["--protein_blastp_mode", "pairwise", "--protein_blastp_max_hits", "3",
         "--protein_blastp_candidate_limit", "4", "--losatp_threads", "2"],
        ["--losat", "losatp", "--losatp_mode", "pairwise", "--losatp_max_hits", "3",
         "--losatp_max_target_seqs", "4", "--losat_threads", "2"],
    ),
    (
        ["--protein_blastp_mode", "orthogroup", "--losatp_bin", "custom-losat"],
        ["--losat", "losatp", "--losatp_mode", "similarity_groups",
         "--losat_bin", "custom-losat"],
    ),
    (
        ["--protein_blastp_mode", "collinear", "--ncbi_blastp_bin", "/opt/blastp",
         "--collinear_min_anchors", "2"],
        ["--losat", "losatp", "--losatp_mode", "collinear",
         "--ncbi_blast_bin", "/opt/blastp", "--collinear_min_anchors", "2"],
    ),
    (["--protein_blastp_mode", "none"], []),
)


def _read_session(path: Path) -> dict:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as handle:
        return json.load(handle)


def _capture_cli_options(monkeypatch, tmp_path, argv):
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    def record(record_id):
        item = SeqRecord(Seq("ATG" * 30), id=record_id, name=record_id)
        item.annotations["molecule_type"] = "DNA"
        item.features.append(
            SeqFeature(FeatureLocation(0, 30, strand=1), type="CDS",
                       qualifiers={"locus_tag": [f"{record_id}_1"]})
        )
        return item

    records = [record("record_a"), record("record_b")]
    captured: dict[str, object] = {}
    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_a, **_k: records)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _p: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _p: None)

    def fake_render(canonical_request, **_kwargs):
        resolved = request_render_module.resolve_request(canonical_request)
        captured["request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
            annotation_warnings=(),
            feature_identity_notices=(),
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render)
    linear_cli_module.linear_main(
        ["--gbk", "a.gb", "b.gb", *argv, "--format", "svg", "-o",
         str(tmp_path / "out"), "--overwrite"]
    )
    request = captured["request"]
    assert isinstance(request, LinearDiagramRequest)
    return request.options


@pytest.mark.parametrize(("argv", "replacement"), RETIRED_FLAGS)
def test_retired_linear_flag_exits_2_and_names_replacement(argv, replacement, capsys):
    with pytest.raises(SystemExit) as excinfo:
        linear_cli_module._get_args(["--gbk", "a.gb", "b.gb", *argv])

    assert excinfo.value.code == 2
    message = capsys.readouterr().err
    assert argv[0] in message
    assert replacement in message


@pytest.mark.parametrize(("legacy", "expected"), LEGACY_REWRITES)
def test_legacy_session_argv_rewrites_retired_flags(legacy, expected):
    from gbdraw.session_io import _canonicalize_legacy_session_cli_args

    args, index_map = _canonicalize_legacy_session_cli_args(
        ["--gbk", "a.gb", *legacy, "-o", "out"],
        mode="linear",
    )

    assert args == ["--gbk", "a.gb", *expected, "-o", "out"]
    # Tokens before and after the rewrite keep their binding positions.
    assert index_map[1] == 1
    assert index_map[2 + len(legacy)] == 2 + len(expected)


def test_release_0_13_0_bgc_session_argv_replays_with_new_flags(tmp_path):
    from gbdraw.session_io import session_to_cli_args

    session = _read_session(FIXTURES / "BGC0000708-BGC0000713.v30.gbdraw-session.json.gz")
    assert session["version"] == 30
    spec = session_to_cli_args(
        session,
        mode="linear",
        temp_dir=tmp_path,
        output_override=None,
        format_override=None,
    )

    for args in (spec.args, spec.cli_invocation_args):
        joined = " ".join(args)
        assert "--losat losatp --losatp_mode similarity_groups" in joined
        assert "--losatp_max_hits 5" in joined
        assert "--similarity_alignment_feature og_1" in joined
        assert "--protein_blastp" not in joined
        assert "--align_orthogroup_feature" not in joined
    parsed = linear_cli_module._get_args(list(spec.args))
    assert parsed.losat == "losatp"
    assert parsed.losatp_mode == "similarity_groups"
    assert parsed.losatp_max_hits == 5
    assert parsed.similarity_alignment_feature == "og_1"


def test_release_0_13_0_protein_session_replays_unchanged_svg(tmp_path, monkeypatch):
    session_path = tmp_path / "legacy.gbdraw-session.json"
    session_path.write_text(
        json.dumps(_read_session(FIXTURES / "cli-linear-protein.v30.gbdraw-session.json.gz")),
        encoding="utf-8",
    )
    monkeypatch.chdir(tmp_path)

    linear_cli_module.linear_main(
        ["--session", str(session_path), "-o", "replay", "-f", "svg", "--save_session"]
    )

    expected = (FIXTURES / "cli-linear-protein.v30.replay.svg").read_bytes()
    assert (tmp_path / "replay.svg").read_bytes() == expected
    saved = json.loads((tmp_path / "replay.gbdraw-session.json").read_text(encoding="utf-8"))
    args = saved["cliInvocation"]["args"]
    assert args[args.index("--losat") : args.index("--losat") + 4] == [
        "--losat", "losatp", "--losatp_mode", "similarity_groups",
    ]
    assert "--protein_blastp_mode" not in args


@pytest.mark.parametrize(("old_argv", "new_argv"), EQUIVALENT_ARGV)
def test_new_flags_map_to_the_typed_options_of_the_old_flags(
    old_argv, new_argv, monkeypatch, tmp_path
):
    from gbdraw.session_io import _canonicalize_legacy_session_cli_args

    rewritten, _ = _canonicalize_legacy_session_cli_args(old_argv, mode="linear")
    assert rewritten == new_argv
    assert _capture_cli_options(monkeypatch, tmp_path, rewritten) == _capture_cli_options(
        monkeypatch, tmp_path, new_argv
    )


def test_new_cli_flags_reach_the_typed_options(monkeypatch, tmp_path):
    from gbdraw.api.options import LosatRuntimeOptions, LosatSearchOptions

    options = _capture_cli_options(
        monkeypatch,
        tmp_path,
        [
            "--losat", "losatp", "--losatp_mode", "collinear",
            "--ncbi_blast_bin", "/opt/ncbi/bin/blastp", "--losat_threads", "6",
            "--losatp_max_hits", "9", "--losatp_max_target_seqs", "123",
            "--losatp_member_max_hits", "3", "--collinear_infer_orthogroups", "off",
        ],
    )

    assert options.losat_search == LosatSearchOptions(
        program="losatp",
        losatp_mode="collinear",
        losatp_max_hits=9,
        losatp_max_target_seqs=123,
        losatp_member_max_hits=3,
        runtime=LosatRuntimeOptions(
            ncbi_blast_executable="/opt/ncbi/bin/blastp",
            threads=6,
        ),
    )
    assert options.collinear_infer_orthogroups is False


def test_cli_defaults_keep_current_losatp_behavior(monkeypatch, tmp_path):
    options = _capture_cli_options(monkeypatch, tmp_path, ["--losat", "losatp"])

    search = options.losat_search
    assert search.losatp_mode == "similarity_groups"
    assert search.losatp_max_hits == 5
    assert search.losatp_max_target_seqs is None
    assert search.losatp_member_max_hits is None
    assert search.runtime.losat_executable is None
    assert options.collinear_infer_orthogroups is True
    assert _capture_cli_options(monkeypatch, tmp_path, []).losat_search is None


@pytest.mark.parametrize(
    "argv",
    (
        ["--losatp_mode", "pairwise"],
        ["--losat", "losatp", "--losat_bin", "a", "--ncbi_blast_bin", "b"],
        ["--losat", "losatp", "--losatp_mode", "pairwise",
         "--similarity_alignment_feature", "x"],
        ["--losat_output_dir", "raw"],
        ["--losat", "losatp", "-b", "a_b.tsv"],
        ["--losat", "losatp", "--losatp_member_max_hits", "0"],
    ),
)
def test_invalid_new_losat_flags_exit_2(argv):
    with pytest.raises(SystemExit) as excinfo:
        linear_cli_module._get_args(["--gbk", "a.gb", "b.gb", *argv])
    assert excinfo.value.code == 2


def test_losat_output_dir_writes_the_raw_losatp_evidence(tmp_path, monkeypatch):
    session = _read_session(FIXTURES / "cli-linear-protein.v30.gbdraw-session.json.gz")
    from gbdraw.session_io import materialize_embedded_file

    paths = [
        str(materialize_embedded_file(item["gb"], temp_dir=tmp_path, role=f"P{index}"))
        for index, item in enumerate(session["files"]["linearSeqs"], start=1)
    ]
    monkeypatch.chdir(tmp_path)
    linear_cli_module.linear_main(
        ["--gbk", *paths, "--losat", "losatp", "-o", "run",
         "--losat_output_dir", "raw", "-f", "svg"]
    )

    raw = (tmp_path / "raw" / "losatp.raw.tsv").read_text(encoding="utf-8")
    assert raw.startswith("# gbdraw raw protein-search evidence\n")
    assert "\n# entry 1: " in raw


def test_typed_options_reject_retired_fields():
    from gbdraw.api.options import LinearDiagramOptions

    for name, value in (
        ("protein_blastp_mode", "pairwise"),
        ("protein_comparison_pairs", ((0, 1),)),
        ("losatp_bin", "losat"),
        ("ncbi_blastp_bin", "blastp"),
        ("losatp_threads", 2),
        ("protein_blastp_max_hits", 5),
        ("protein_blastp_candidate_limit", 4),
        ("orthogroup_member_max_hits", 3),
    ):
        with pytest.raises(TypeError, match=name):
            LinearDiagramOptions(**{name: value})


def test_losat_search_options_validate_program_and_mode():
    from gbdraw.api.options import LosatRuntimeOptions, LosatSearchOptions
    from gbdraw.exceptions import ValidationError

    assert LosatSearchOptions(program="losatn").losatn_task == "megablast"
    with pytest.raises(ValidationError, match="losatp_mode"):
        LosatSearchOptions(program="losatp")
    with pytest.raises(ValidationError, match="pairs"):
        LosatSearchOptions(program="losatp", losatp_mode="collinear", pairs=((0, 1),))
    with pytest.raises(ValidationError, match="losatp_max_hits"):
        LosatSearchOptions(program="losatp", losatp_mode="pairwise", losatp_max_hits=0)
    with pytest.raises(ValidationError):
        LosatRuntimeOptions(losat_executable="a", ncbi_blast_executable="b")
    assert LosatRuntimeOptions(losat_executable="losat").losat_executable is None
    search = LosatSearchOptions(
        program="losatp", losatp_mode="pairwise", pairs=[[0, 1]]
    )
    assert search.pairs == ((0, 1),)


def test_introductory_api_uses_the_new_comparison_fields():
    from gbdraw import LinearComparisonOptions, LinearOptions
    from gbdraw.api.options import LosatRuntimeOptions, LosatSearchOptions
    from gbdraw.interface import _linear_options

    defaults = LinearComparisonOptions()
    assert defaults.losat is None
    assert defaults.losatp_mode == "similarity_groups"
    assert defaults.losat_executable is None
    assert defaults.ncbi_blast_executable is None
    assert defaults.max_target_seqs is None
    assert defaults.member_max_hits is None
    assert _linear_options(LinearOptions(), record_count=2).losat_search is None

    options = _linear_options(
        LinearOptions(
            comparisons=LinearComparisonOptions(
                losat="losatp",
                losatp_mode="pairwise",
                pairs=[(0, 1)],
                max_hits=4,
                max_target_seqs=8,
                member_max_hits=2,
                ncbi_blast_executable="/opt/blastp",
                threads=3,
            )
        ),
        record_count=2,
    )
    assert options.losat_search == LosatSearchOptions(
        program="losatp",
        losatp_mode="pairwise",
        pairs=((0, 1),),
        losatp_max_hits=4,
        losatp_max_target_seqs=8,
        losatp_member_max_hits=2,
        runtime=LosatRuntimeOptions(ncbi_blast_executable="/opt/blastp", threads=3),
    )


@pytest.mark.parametrize(
    "name",
    ("protein_mode", "blastp_executable", "candidate_limit", "orthogroup_member_max_hits"),
)
def test_introductory_api_rejects_retired_fields(name):
    from gbdraw import LinearComparisonOptions

    with pytest.raises(TypeError, match=name):
        LinearComparisonOptions(**{name: None})


def _protein_marker(payload: dict) -> dict:
    markers = [
        item for item in payload["comparisons"]
        if item["kind"] == "generatedProteinComparison"
    ]
    assert len(markers) == 1
    return markers[0]


def _round_trip(session_path: Path, tmp_path: Path):
    from gbdraw.session import load_session_document, materialize_session, session_to_request
    from gbdraw.session_request_codec import encode_canonical_request

    session = _read_session(session_path)
    document = load_session_document(session)
    with materialize_session(document, output_directory=tmp_path) as materialized:
        request = session_to_request(materialized)
        encoded = encode_canonical_request(request)
    return session, request, encoded.payload


@pytest.mark.parametrize(
    "fixture",
    (
        "q-frame-main-linear-reverse.v42.gbdraw-session.json.gz",
        "se06-main-linear-blast-cli.v42.gbdraw-session.json.gz",
        "se08-main-linear-cli.v42.gbdraw-session.json.gz",
    ),
)
def test_schema_7_protein_marker_round_trips_wire_names(fixture, tmp_path):
    session, request, payload = _round_trip(FIXTURES / fixture, tmp_path)

    assert session["renderRequest"]["schema"] == 7
    assert request.options.losat_search is None
    before = copy.deepcopy(_protein_marker(session["renderRequest"]))
    # Schemas 1-7 also stored the retired alignment target; schema 8 does not.
    assert before["settings"].pop("alignOrthogroupFeature") is None
    assert _protein_marker(payload) == before


@pytest.mark.parametrize(
    ("session_path", "mode", "threads"),
    (
        # The Gallery Session stores its collinear result (mode "none").
        (GALLERY_SESSIONS / "vibrio-harveyi-group-collinear.gbdraw-session.json.gz",
         "none", 16),
        (GALLERY_SESSIONS / "hepatoplasmataceae_orthogroup.gbdraw-session.json.gz",
         "none", 32),
    ),
)
def test_protein_marker_settings_round_trip_through_typed_options(
    session_path, mode, threads, tmp_path
):
    session, request, payload = _round_trip(session_path, tmp_path)

    search = request.options.losat_search
    assert search.program == "losatp"
    assert search.losatp_mode == mode
    assert search.runtime.threads == threads
    before = copy.deepcopy(_protein_marker(session["renderRequest"]))
    assert _protein_marker(payload) == before
