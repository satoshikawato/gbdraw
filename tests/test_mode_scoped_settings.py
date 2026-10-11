"""Session 46 mode slices: the registry, its Web leaf, the split twin, and the 46 validators."""

from __future__ import annotations

import copy
import csv
import gzip
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

import pytest

import gbdraw.cli_utils.session as cli_session_module
from gbdraw.api.session_compat import _session_protein_mode
from gbdraw.exceptions import ValidationError
from gbdraw.linear import linear_main
from gbdraw.session import load_session_document
from gbdraw.session_io import (
    CURRENT_SESSION_VERSION,
    FLAT_DRAFT_TOP_LEVEL_FIELDS,
    SUPPORTED_SESSION_VERSIONS,
    expand_session_feature_catalog,
    migrate_session_flat_draft,
    session_depth_source_widths,
    mode_split_palette_colors,
    session_mode,
    split_draft_into_modes,
    validate_session,
)
from gbdraw.web_support.mode_scoped_settings import (
    DIAGRAM_MODES,
    MIGRATION_TOKENS,
    MODE_SCOPED_SETTINGS,
    MODE_SCOPED_SETTINGS_REVISION,
    SLICE_CONTAINERS,
    unmanaged_config_override_modes,
)

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]
FIXTURES = Path(__file__).parent / "fixtures" / "sessions"
FIELDS_INVALID = {"code": "INPUT_INVALID", "field": "schema", "reason": "FIELDS"}


def _read(name: str) -> dict[str, Any]:
    return json.loads(gzip.decompress((FIXTURES / name).read_bytes()))


def _upgrade(session: dict[str, Any]) -> dict[str, Any]:
    """A Session 27-44 as the CLI writes its draft: the older migrations, then the split."""

    session = expand_session_feature_catalog(session)
    migrated = migrate_session_flat_draft(session).session
    return split_draft_into_modes(migrated, committed_mode=session_mode(session))


# --- registry and the Web leaf ------------------------------------------------


def test_registry_rows_are_unique_and_use_known_tokens() -> None:
    names = [(row.domain, row.path) for row in MODE_SCOPED_SETTINGS]
    assert len(names) == len(set(names))
    assert {row.migrate for row in MODE_SCOPED_SETTINGS} <= MIGRATION_TOKENS
    for row in MODE_SCOPED_SETTINGS:
        assert row.modes in ("both", *DIAGRAM_MODES), row
        # A row that goes to its own mode names that mode.
        assert (row.migrate == "own") == (row.modes != "both" and row.migrate != "profile"), row
    assert {row.key for row in MODE_SCOPED_SETTINGS if row.key} == {
        "id", "JSON[sourceUid,selector]", "JSON[recordKey,biologicalFeatureId]", "recordKey\\0featureId",
    }
    assert SLICE_CONTAINERS[""] == {"config", "features", "editorState", "ui"}
    assert SLICE_CONTAINERS["config.losat"] == {"outfmt", "blastn", "blastp"}


def test_registry_rows_match_the_phase_e_registry_table() -> None:
    # The table (Phase E's registry-v46.tsv, revision 5) is the readable form
    # of the registry: each of its slice rows is a row here, and the reverse.
    assert MODE_SCOPED_SETTINGS_REVISION == 5
    with (FIXTURES / "mode-scoped-settings-registry-v46.tsv").open(encoding="utf-8", newline="") as handle:
        table = list(csv.DictReader(handle, delimiter="\t"))
    assert len(table) == 278
    rows = set()
    for row in table:
        if row["migrate"] == "-":
            continue
        path, modes = row["registry_entry"].split(" | modes=")
        assert modes == row["modes"], row["item"]
        # One table row names both pending palette keys.
        paths = (
            ("ui.pendingPaletteName", "ui.pendingPaletteColors")
            if path == "ui.pendingPaletteName/Colors"
            else (path,)
        )
        rows |= {(name, row["modes"], row["key"] or None, row["migrate"]) for name in paths}

    assert rows == {(f"{row.domain}.{row.path}", row.modes, row.key, row.migrate) for row in MODE_SCOPED_SETTINGS}


def test_generated_web_mode_scoped_settings_match_python_source() -> None:
    result = subprocess.run(
        [sys.executable, "tools/generate_mode_scoped_settings.py", "--check"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


def test_unmanaged_override_leaves_belong_to_the_modes_that_render_them() -> None:
    assert unmanaged_config_override_modes("objects.ticks.tick_width") == ("circular",)
    assert unmanaged_config_override_modes("objects.blast_match.curve_tension") == ("linear",)
    assert unmanaged_config_override_modes("canvas.dpi") == ("circular", "linear")


# --- the split ----------------------------------------------------------------


def _draft(**config: Any) -> dict[str, Any]:
    return {"ui": {"mode": "circular"}, "config": config}


def _split(draft: dict[str, Any], **context: Any) -> dict[str, Any]:
    context.setdefault("committed_mode", "circular")
    context.setdefault("depth_sources", {"circular": 0, "linear": 0})
    return split_draft_into_modes(draft, **context)


def test_split_copies_shared_values_and_gives_mode_only_values_to_their_mode() -> None:
    result = _split(_draft(
        form={"prefix": "p", "species": "S", "show_gc": True},
        adv={"nt": "AT", "feature_height": 30},
        palette="dark",
    ))

    circular, linear = (result["modes"][mode]["config"] for mode in DIAGRAM_MODES)
    assert circular["form"] == {"prefix": "p", "species": "S"}
    assert linear["form"] == {"prefix": "p", "show_gc": True}
    # The shown mode also keeps the ribbons a draft without a match style drew.
    assert circular["adv"] == {"nt": "AT", "pairwise_match_style": "ribbon"}
    assert linear["adv"] == {"nt": "AT", "feature_height": 30}
    assert circular["palette"] == linear["palette"] == "dark"
    assert set(FLAT_DRAFT_TOP_LEVEL_FIELDS).isdisjoint(result)


def test_split_gives_the_shown_mode_its_flat_profile_values_and_the_other_its_profile() -> None:
    profiles = {
        "schema": 1,
        "activeMode": "linear",
        "profiles": {
            "circular": {"values": {"evalue": "1e-3", "plot_title": "C"}, "managed": {}},
            "linear": {"values": {"evalue": "stale", "plot_title": "stale"}, "managed": {}},
        },
    }
    draft = {
        "ui": {"mode": "linear"},
        "config": {"form": {"plot_title": "L"}, "adv": {"evalue": "1e-2"}, "modeProfiles": profiles},
    }

    result = _split(draft, committed_mode="circular")

    circular, linear = (result["modes"][mode]["config"] for mode in DIAGRAM_MODES)
    assert (linear["form"]["plot_title"], linear["adv"]["evalue"]) == ("L", "1e-2")
    assert (circular["form"]["plot_title"], circular["adv"]["evalue"]) == ("C", "1e-3")
    # No saved value: that mode's default (absent). A draft without a pairwise
    # match style drew ribbons in its shown mode.
    assert "identity" not in circular["adv"]
    assert linear["adv"]["pairwise_match_style"] == "ribbon"
    assert "modeProfiles" not in circular and "modeProfiles" not in linear


def test_split_moves_layout_slots_depth_and_show_depth_per_mode() -> None:
    draft = {
        "ui": {
            "mode": "circular",
            "layoutPreferences": {
                "circular": {"single": {"legend": "left"}, "multi": {"legend": "right"}},
                "linear": {"legend": "top", "plotTitlePosition": "bottom"},
            },
        },
        "config": {
            "form": {"show_depth": True},
            "adv": {
                "depth_large_tick_interval": 9,
                "depth_tracks": [{"label": "a", "large_tick_interval": None}, {"label": "b"}],
            },
        },
    }

    result = _split(draft, depth_sources={"circular": 2, "linear": 0})

    circular, linear = (result["modes"][mode] for mode in DIAGRAM_MODES)
    assert circular["ui"]["layoutPreferences"] == {"single": {"legend": "left"}, "multi": {"legend": "right"}}
    assert linear["ui"]["layoutPreferences"] == {"legend": "top", "plotTitlePosition": "bottom"}
    assert "layoutPreferences" not in result["ui"]
    assert circular["config"]["form"]["show_depth"] is True
    assert linear["config"]["form"]["show_depth"] is False
    # Each slice keeps its own sources' series, and the flat tick value that a
    # series without its own reads.
    assert circular["config"]["adv"]["depth_tracks"] == [
        {"label": "a", "large_tick_interval": None},
        {"label": "b"},
    ]
    assert linear["config"]["adv"]["depth_tracks"] == [{"label": "a", "large_tick_interval": None}]
    assert circular["config"]["adv"]["depth_large_tick_interval"] == 9
    assert linear["config"]["adv"]["depth_large_tick_interval"] == 9


def test_split_gives_the_saved_results_mode_the_legend_and_feature_edits() -> None:
    triple = json.dumps(["linear", "record-1", "f1"], separators=(",", ":"))
    row = {"recordKey": "record-1", "biologicalFeatureId": "f1", "featureVisibility": "off",
           "labelVisibility": None, "labelText": None, "labelSourceText": None}
    draft = {
        "ui": {"mode": "circular", "canvasPadding": {"top": 1}},
        "config": {},
        "features": {"featureOverrides": {triple: {"scope": "linear", **row}}, "labelTextBulkOverrides": {"a": "b"},
                     "selectedFeatureRecordIdx": 0},
        "editorState": {
            "legend": {"entries": [{"caption": "CDS"}], "originalOrder": ["CDS"]},
            "featureStrokes": {"overrides": {}},
            "originalSvgStroke": {"color": "gray", "width": 1},
        },
    }

    result = _split(draft, committed_mode="linear")

    circular, linear = (result["modes"][mode] for mode in DIAGRAM_MODES)
    pair = json.dumps(["record-1", "f1"], separators=(",", ":"))
    assert linear["features"]["featureOverrides"] == {pair: row}
    assert "featureOverrides" not in circular["features"]
    assert circular["features"]["labelTextBulkOverrides"] == {"a": "b"}
    assert linear["editorState"] == {"legend": {"entries": [{"caption": "CDS"}]}, "featureStrokes": {"overrides": {}}}
    assert "editorState" not in circular
    # The original SVG stroke is read as Web Load reads it: a color name as hex (OV-160).
    assert result["editorState"] == {"legend": {"originalOrder": ["CDS"]}, "originalSvgStroke": {"color": "#808080", "width": 1}}
    assert circular["ui"]["canvasPadding"] == linear["ui"]["canvasPadding"] == {"top": 1}
    assert "canvasPadding" not in result["ui"]
    assert "features" not in result


def test_split_places_scoped_rows_and_override_leaves_in_their_modes() -> None:
    def placement(scope: str, side: str | None) -> tuple[str, dict[str, Any]]:
        target = {"kind": "main"} if side is None else {"kind": "lane", "side": side, "level": 1}
        row = {"scope": scope, "recordKey": "r", "biologicalFeatureId": f"f-{side}", "placement": target}
        return json.dumps([scope, "r", f"f-{side}"], separators=(",", ":")), row

    draft = _draft(
        recordDisplayDrafts=[{"scope": "linear", "sourceUid": "u", "selector": "#1"}],
        featurePlacementOverrides=dict(
            placement(*args) for args in (("circular", "outward"), ("circular", None), ("linear", None))
        ),
        unmanagedConfigOverrides={"objects.ticks.tick_width": 4, "objects.blast_match.curve_tension": 0.3,
                                  "canvas.dpi": 300},
    )

    result = _split(draft)

    circular, linear = (result["modes"][mode]["config"] for mode in DIAGRAM_MODES)
    assert circular["recordDisplayDrafts"] == []
    assert linear["recordDisplayDrafts"] == [{"sourceUid": "u", "selector": "#1"}]
    assert set(circular["featurePlacementOverrides"]) == {
        json.dumps(["r", "f-outward"], separators=(",", ":")), json.dumps(["r", "f-None"], separators=(",", ":")),
    }
    assert list(linear["featurePlacementOverrides"]) == [json.dumps(["r", "f-None"], separators=(",", ":"))]
    assert all("scope" not in row for row in linear["featurePlacementOverrides"].values())
    assert circular["unmanagedConfigOverrides"] == {"objects.ticks.tick_width": 4, "canvas.dpi": 300}
    assert linear["unmanagedConfigOverrides"] == {"objects.blast_match.curve_tension": 0.3, "canvas.dpi": 300}


def test_split_keeps_an_annotation_bound_to_one_mode_in_that_mode() -> None:
    linear_binding = '["linear","seq-1","gb",["a.gb",1],null]::[0,"A",10]'
    annotations = [
        {"id": "bound", "target": {"kind": "featureSpan"}, "metadata": {"_gbdraw_web_target_record_key": linear_binding}},
        {"id": "selected", "target": {"kind": "featureIdentity", "scope": "circular", "recordKey": "r",
                                      "biologicalFeatureId": "f"}, "metadata": {}},
        {"id": "free", "target": {"kind": "coordinateSpan"}, "metadata": {}},
    ]

    result = _split(_draft(annotationSets=[{"id": "set", "annotations": annotations}]))

    circular, linear = (result["modes"][mode]["config"]["annotationSets"] for mode in DIAGRAM_MODES)
    assert [item["id"] for item in circular[0]["annotations"]] == ["selected", "free"]
    assert [item["id"] for item in linear[0]["annotations"]] == ["bound", "free"]
    assert "scope" not in circular[0]["annotations"][0]["target"]


def test_split_moves_losat_execution_to_the_app_and_blastp_to_linear() -> None:
    losat = {"outfmt": "6", "executionMode": "threaded", "threadsPerJob": "auto", "blastn": {"task": "megablast"},
             "blastp": {"mode": "collinear", "candidateLimit": 7,
                        "hitLimitsByMode": {"collinear": {"candidateLimit": 5}}}}

    result = _split(_draft(
        losat=losat, paletteInstantPreviewEnabled=True, adv={"rich_feature_popup": False},
        cliOptions={"rawArgs": ["--gbk", "a.gb"]},
    ))

    # App-level settings and Session provenance leave the slices.
    assert result["ui"]["richFeaturePopup"] is False
    assert result["cliOptions"] == {"rawArgs": ["--gbk", "a.gb"]}
    for mode in DIAGRAM_MODES:
        assert "rich_feature_popup" not in result["modes"][mode]["config"].get("adv", {})
        assert "cliOptions" not in result["modes"][mode]["config"]
    circular, linear = (result["modes"][mode]["config"]["losat"] for mode in DIAGRAM_MODES)
    assert circular == {"outfmt": "6", "blastn": {"task": "megablast"}}
    # PD-OI-002: blastp moves unchanged; the flat limit stays the effective one.
    assert linear["blastp"] == losat["blastp"]
    assert result["ui"]["losatExecution"] == {"executionMode": "threaded", "threadsPerJob": "auto"}
    assert result["ui"]["paletteInstantPreviewEnabled"] is True


def test_split_counts_depth_sources_in_the_web_file_bindings() -> None:
    binding = {"resourceId": "d"}
    assert session_depth_source_widths(
        {"c_depth": [[binding, None]], "linearSeqs": [{"depth": [None]}, {"depth": [None, binding]}]}
    ) == {"circular": 2, "linear": 2}
    assert session_depth_source_widths({"c_depth": None, "linearSeqs": [{"depth": [None]}]}) == {
        "circular": 0, "linear": 0,
    }


def test_split_needs_the_saved_results_mode_and_writes_no_slices_without_a_draft() -> None:
    with pytest.raises(ValidationError, match="saved Result's mode") as excinfo:
        split_draft_into_modes({"config": {}}, committed_mode=None)
    assert excinfo.value.diagnostic == FIELDS_INVALID
    assert "modes" not in _split({"ui": {"mode": "linear", "canvasPadding": {}}})


# --- Session 46 validation ----------------------------------------------------


@pytest.fixture(scope="module")
def written_session_46(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """A Session 46 that the CLI wrote from a Web-saved two-mode Session 44."""

    work = tmp_path_factory.mktemp("session-46")
    sidecar = work / "replay.gbdraw-session.json"
    linear_main(
        ["--session", str(FIXTURES / "feature-placements-linear.v44.gbdraw-session.json.gz"),
         "--output", str(work / "replay"), "--format", "svg", "--session_output", str(sidecar)]
    )
    return load_session_document(sidecar).to_dict()


@pytest.fixture
def session_46(written_session_46: dict[str, Any]) -> dict[str, Any]:
    return copy.deepcopy(written_session_46)


def test_session_46_reads_two_mode_slices_and_rejects_session_45(session_46: dict[str, Any]) -> None:
    assert CURRENT_SESSION_VERSION == 46 and 45 not in SUPPORTED_SESSION_VERSIONS
    assert set(session_46["modes"]) == {"circular", "linear"}
    validate_session(session_46)
    with pytest.raises(ValidationError, match="Unsupported session version: 45"):
        validate_session({**session_46, "version": 45})
    with pytest.raises(ValidationError, match="Session version 44 cannot contain modes"):
        validate_session({**session_46, "version": 44})


@pytest.mark.parametrize(
    ("damage", "message"),
    [
        (lambda s: s.update(config={}), "it cannot contain config"),
        (lambda s: s.update(features={}), "it cannot contain features"),
        (lambda s: s["editorState"].setdefault("legend", {}).update(entries=[]),
         "it cannot contain editorState.legend.entries"),
        (lambda s: s["editorState"].update(featureStrokes={"overrides": {}}),
         "it cannot contain editorState.featureStrokes"),
        (lambda s: s["ui"].update(layoutPreferences={}), "it cannot contain ui.layoutPreferences"),
        (lambda s: s["ui"].update(losatExecution={"modeProfiles": 1}), "ui.losatExecution holds only"),
        (lambda s: s["modes"].update(third={}), "a circular and a linear slice only"),
        (lambda s: s["modes"]["circular"].update(results=[]), "modes.circular cannot contain results"),
        (lambda s: s["modes"]["linear"]["config"].update(modeProfiles={}),
         "modes.linear.config cannot contain modeProfiles"),
        (lambda s: s["modes"]["linear"]["config"]["losat"].update(executionMode="auto"),
         "modes.linear.config.losat cannot contain executionMode"),
        (lambda s: s["modes"]["linear"]["config"]["form"].update(legend="top"),
         "modes.linear.config.form cannot contain legend"),
    ],
)
def test_session_46_rejects_draft_fields_outside_the_mode_slices(
    session_46: dict[str, Any], damage: Any, message: str
) -> None:
    damaged = copy.deepcopy(session_46)
    damage(damaged)

    with pytest.raises(ValidationError, match=message) as excinfo:
        validate_session(damaged)
    assert excinfo.value.diagnostic == FIELDS_INVALID


@pytest.mark.parametrize(
    ("mode", "path", "value"),
    [("circular", "objects.ticks.tick_width", 4), ("linear", "objects.blast_match.curve_tension", 0.3)],
)
def test_session_46_validates_unmanaged_overrides_against_their_slice_mode(
    session_46: dict[str, Any], mode: str, path: str, value: Any
) -> None:
    # OV-106: each mode keeps the overrides of its own renderer.
    session = session_46
    session["modes"][mode]["config"]["unmanagedConfigOverrides"] = {path: value}
    validate_session(session)

    other = "linear" if mode == "circular" else "circular"
    session["modes"][other]["config"]["unmanagedConfigOverrides"] = {path: value}
    with pytest.raises(ValidationError, match="cannot target") as excinfo:
        validate_session(session)
    assert excinfo.value.diagnostic is not None
    assert excinfo.value.diagnostic["code"] == "MODE_SETTING"


def test_session_46_settings_only_needs_the_shown_modes_configuration() -> None:
    source = _read("settings-only.v42.json.gz")
    session = split_draft_into_modes(source, committed_mode=source["ui"]["mode"]) | {"version": 46}
    validate_session(session)
    assert load_session_document(session).has_canonical_request is False

    shown = session["ui"]["mode"]
    session["modes"][shown]["config"].pop("form")
    with pytest.raises(ValidationError, match="active Web configuration"):
        validate_session(session)


def test_session_protein_mode_reads_the_linear_slice() -> None:
    blastp = {"config": {"losat": {"blastp": {"mode": "collinear"}}}}
    assert _session_protein_mode({"modes": {"linear": blastp}}) == "collinear"
    assert _session_protein_mode({"modes": {"circular": blastp}}) is None
    assert _session_protein_mode(blastp) == "collinear"


# --- the CLI re-save ----------------------------------------------------------


def test_cli_resave_keeps_a_session_46_slices(tmp_path: Path) -> None:
    fixture = FIXTURES / "feature-placements-linear.v44.gbdraw-session.json.gz"
    first = tmp_path / "first.gbdraw-session.json"
    second = tmp_path / "second.gbdraw-session.json"
    for source, sidecar in ((fixture, first), (first, second)):
        linear_main(
            ["--session", str(source), "--output", str(tmp_path / sidecar.stem), "--format", "svg",
             "--session_output", str(sidecar)]
        )

    written = [load_session_document(path).to_dict() for path in (first, second)]
    assert [payload["version"] for payload in written] == [46, 46]
    assert written[1]["modes"] == written[0]["modes"]
    assert "config" not in written[1] and "features" not in written[1]


def test_cli_resave_of_a_config_less_session_44_moves_its_request_annotation_targets(
    tmp_path: Path, caplog: pytest.LogCaptureFixture
) -> None:
    # OV-135: a Session 40-44 written without a Web draft (CLI, Python API)
    # takes its request's annotation sets as its draft, as Web Load does, so
    # its certain hash= target moves to the feature's source identity.
    source = _read("selected-feature-annotations.v44.gbdraw-session.json.gz")
    source.pop("config")
    source_path = tmp_path / "config-less.v44.gbdraw-session.json"
    source_path.write_text(json.dumps(source), encoding="utf-8")
    sidecar = tmp_path / "replay.gbdraw-session.json"
    caplog.set_level("INFO", logger=cli_session_module.__name__)

    linear_main(
        ["--session", str(source_path), "--output", str(tmp_path / "replay"), "--format", "svg",
         "--session_output", str(sidecar)]
    )

    saved = load_session_document(sidecar).to_dict()
    targets = [
        (annotation["id"], annotation["target"]["kind"])
        for annotation_set in saved["modes"]["linear"]["config"]["annotationSets"]
        for annotation in annotation_set["annotations"]
    ]
    assert [kind for _, kind in targets] == ["featureIdentity", "featureSpan", "featureSpan"]
    assert set(saved["modes"]["linear"]["config"]) == {"annotationSets"}
    # Every target is bound to a Linear record.
    assert all(not annotation_set["annotations"] for annotation_set in saved["modes"]["circular"]["config"]["annotationSets"])
    assert saved["renderRequest"]["diagramOptions"]["annotations"] == source["renderRequest"]["diagramOptions"]["annotations"]
    assert "INFO: 1 annotation(s) from Session version 44 named a feature by hash=" in caplog.text


def test_split_moves_a_boolean_instant_preview_to_the_app() -> None:
    # Registry row 174, for any source version. Web Load applies ui, then a
    # boolean draft value, so the draft value wins over a saved ui value.
    assert _split(_draft(paletteInstantPreviewEnabled=True))["ui"]["paletteInstantPreviewEnabled"] is True
    shown = {"ui": {"mode": "circular", "paletteInstantPreviewEnabled": False},
             "config": {"paletteInstantPreviewEnabled": True}}
    result = _split(shown)
    assert result["ui"]["paletteInstantPreviewEnabled"] is True
    unread = {"ui": {"mode": "circular", "paletteInstantPreviewEnabled": False},
              "config": {"paletteInstantPreviewEnabled": "yes"}}
    assert _split(unread)["ui"]["paletteInstantPreviewEnabled"] is False
    assert "paletteInstantPreviewEnabled" not in result["modes"]["circular"].get("config", {})


def test_split_resolves_override_colors_against_the_draft_palette() -> None:
    # Web Load reads colors marked colorsAreOverrides over the palette's colors
    # (services/config.js applyConfigData); the split stores the result.
    with (REPO_ROOT / "gbdraw" / "data" / "color_palettes.toml").open("rb") as handle:
        palette = {key: str(value) for key, value in tomllib.load(handle)["default"].items()}
    draft = _draft(palette="default", colors={"CDS": " #123456 ", "collinear_block_2": "#000000"},
                   colorsAreOverrides=True)
    result = _split(draft, palette_colors=mode_split_palette_colors(draft["config"]))

    colors = result["modes"]["circular"]["config"]["colors"]
    assert result["modes"]["linear"]["config"]["colors"] == colors
    assert colors["CDS"] == "#123456"
    assert colors["tRNA"] == palette["tRNA"]
    assert "collinear_block_2" not in colors
    assert colors["pairwise_match"] == palette.get("pairwise_match", "#d3d3d3")
    for mode in DIAGRAM_MODES:
        assert "colorsAreOverrides" not in result["modes"][mode]["config"]

    # Without the flag the stored colors are complete and stay as saved.
    complete = _split(_draft(palette="default", colors={"CDS": "red"}, colorsAreOverrides=False))
    assert complete["modes"]["circular"]["config"]["colors"] == {"CDS": "red"}
    assert "colorsAreOverrides" not in complete["modes"]["circular"]["config"]
