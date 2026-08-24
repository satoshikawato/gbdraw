from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path

import pytest

from tools.capture_gallery_tutorial_screenshots import (
    DEFAULT_MAX_FILE_SIZE_KB,
    DEFAULT_MAX_IMAGE_HEIGHT,
    DEFAULT_MAX_IMAGE_WIDTH,
    TutorialContext,
    ValidationResult,
    add_capture_contract_validation,
    iter_operation_contexts,
    load_ready_examples,
    resolve_gallery_reference,
    validate_tutorial_media,
    wait_for_app_shell,
)


pytestmark = pytest.mark.gallery


BASIC_EXAMPLES = {
    "HmmtDNA_basic_circular": {
        "operation_count": 6,
        "data_dependent_count": 5,
        "filename": "HmmtDNA.gbk",
        "resource_marker": b"VERSION     NC_012920.1",
    },
    "lambda_basic_linear": {
        "operation_count": 9,
        "data_dependent_count": 8,
        "filename": "NC_001416.gb",
        "resource_marker": b"VERSION     NC_001416.1",
    },
}


LAYOUT_RESULT_CAPTURES = {
    "HmmtDNA_basic_circular": {
        "src": "./media/HmmtDNA_basic_circular/manual-04-01-final-preview.webp",
        "selector": "svg[data-capture-hmmt-basic-final='true']",
        "viewport": {"width": 1400, "height": 1000},
        "state": {
            "mode": "circular",
            "form.legend": "right",
            "form.track_type": "middle",
            "form.labels_mode": "out",
        },
        "visible_text": {"Homo sapiens", "NC_012920.1", "GC content", "GC skew (+)"},
    },
    "lambda_basic_linear": {
        "src": "./media/lambda_basic_linear/manual-05-01-final-preview.webp",
        "selector": "svg[data-capture-lambda-final='true']",
        "viewport": {"width": 1900, "height": 900},
        "state": {
            "mode": "linear",
            "linearSeqs.length": 1,
            "linearComparisonPlan.mode": "none",
            "form.legend": "left",
            "form.show_labels_linear": "all",
            "form.scale_style": "ruler",
        },
        "visible_text": {"NC_001416.1", "48,502 bp", "major capsid protein", "CDS"},
    },
    "tobacco-chloroplast": {
        "src": "./media/tobacco-chloroplast/manual-08-01-chloroplast-preview.webp",
        "selector": "svg[data-capture-chloroplast-final='true']",
        "viewport": {"width": 1400, "height": 1000},
        "state": {
            "mode": "circular",
            "form.legend": "upper_left",
            "annotationSets.0.id": "plastome_regions",
            "adv.circular_track_slots_enabled": True,
        },
        "visible_text": {"Nicotiana tabacum", "NC_001879.2", "LSC", "IRb", "SSC", "IRa", "GC content"},
    },
    "Vnig_TUMSAT-TG-2018": {
        "src": "./media/Vnig_TUMSAT-TG-2018/manual-06-01-multirecord-preview.webp",
        "selector": "svg[data-capture-vnig-multirecord-final='true']",
        "viewport": {"width": 1600, "height": 1100},
        "state": {
            "mode": "circular",
            "form.multi_record_canvas": True,
            "form.legend": "left",
            "adv.multi_record_size_mode": "auto",
            "adv.multi_record_min_radius_ratio": 0.55,
            "adv.multi_record_column_gap_ratio": 0.1,
            "adv.multi_record_row_gap_ratio": 0.05,
            "adv.multi_record_positions": [
                {"selector": "#1", "row": 1},
                {"selector": "#2", "row": 1},
                {"selector": "#3", "row": 2},
                {"selector": "#4", "row": 2},
                {"selector": "#5", "row": 2},
                {"selector": "#6", "row": 2},
            ],
            "form.plot_title": "<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018, complete genome",
            "adv.plot_title_position": "bottom",
        },
        "visible_text": {"Vibrio nigripulchritudo", "AP024087.1", "AP024088.1", "pVNTG1", "pVNTG4", "GC content", "GC skew (+)"},
    },
    "vibrio-harveyi-group-collinear": {
        "src": "./media/vibrio-harveyi-group-collinear/manual-08-01-collinear-overview.webp",
        "selector": "svg[data-capture-vibrio-overview='true']",
        "viewport": {"width": 1800, "height": 1100},
        "state": {
            "mode": "linear",
            "linearSeqs.length": 11,
            "linearRecordLayoutEnabled": True,
            "linearRecordGap": 48,
            "linearRecordRows.length": 11,
            "form.legend": "bottom",
            "form.plot_title": "",
            "form.linear_track_layout": "above",
            "form.scale_style": "ruler",
            "form.show_labels_linear": "none",
            "losatProgram": "blastp",
            "losat.blastp.mode": "collinear",
            "losat.blastp.collinearColorMode": "orientation_identity",
        },
        "visible_text": {"Vibrio harveyi", "Vibrio owensii", "Vibrio campbellii", "Vibrio parahaemolyticus", "Vibrio alginolyticus", "Collinear", "Inverted"},
    },
    "HmmtDNA_ATskew": {
        "src": "./media/HmmtDNA_ATskew/manual-09-01-atskew-preview.webp",
        "selector": "svg[data-capture-hmmt-atskew-final='true']",
        "viewport": {"width": 1400, "height": 1000},
        "state": {
            "mode": "circular",
            "form.legend": "left",
            "form.track_type": "middle",
            "form.labels_mode": "out",
            "adv.circular_track_slots_enabled": True,
            "adv.circular_track_slots.length": 5,
            "adv.circular_track_slots.0.id": "features",
            "adv.circular_track_slots.1.id": "gc_content",
            "adv.circular_track_slots.2.id": "gc_skew",
            "adv.circular_track_slots.3.id": "a_skew_2",
            "adv.circular_track_slots.3.params.nt": "AT",
            "adv.circular_track_slots.3.params.legend_label": "AT skew",
            "adv.circular_track_slots.3.params.positive_color": "#deaf6e",
            "adv.circular_track_slots.3.params.negative_color": "#7294e3",
            "adv.circular_track_slots.4.id": "ticks",
            "adv.circular_track_slots.4.params.tick_label_layout": "label_in_tick_out",
        },
        "visible_text": {"Homo sapiens", "NC_012920.1", "GC content", "GC skew (+)", "AT skew (+)", "AT skew (-)"},
    },
}


LINEAR_COMPARISON_PANEL_EXAMPLES = {
    "BGC0000708-BGC0000713",
    "hepatoplasmataceae_orthogroup",
    "hepatoplasmataceae_collinear",
    "majanivirus_orthogroup",
    "vibrio-harveyi-group-collinear",
}


LOSATP_MODE_CAPTURES = {
    "BGC0000708-BGC0000713": {
        "src": "./media/BGC0000708-BGC0000713/manual-03-02-select-losatp-orthogroups.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "hepatoplasmataceae_orthogroup": {
        "src": "./media/hepatoplasmataceae_orthogroup/manual-04-01-orthogroups-mode.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "hepatoplasmataceae_collinear": {
        "src": "./media/hepatoplasmataceae_collinear/manual-04-01-collinear-reduction.webp",
        "selected": "Collinear blocks",
        "state": "collinear",
    },
    "majanivirus_orthogroup": {
        "src": "./media/majanivirus_orthogroup/manual-03-02-losatp-orthogroups.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "vibrio-harveyi-group-collinear": {
        "src": "./media/vibrio-harveyi-group-collinear/manual-04-01-search-method-collinear.webp",
        "selected": "Collinear blocks",
        "state": "collinear",
    },
}


FIRST_COMPARISON_BOUNDARY_CAPTURE = {
    "src": "./media/BGC0000708-BGC0000713/manual-03-03-first-comparison-boundary.webp",
    "selector": "[data-edge-key='record-1->record-2']",
    "viewport": {"width": 1600, "height": 1100},
    "state": {
        "mode": "linear",
        "linearSeqs.0.uid": "record-1",
        "linearSeqs.1.uid": "record-2",
        "linearComparisonResolution.edges.0.edgeKey": "record-1->record-2",
        "linearComparisonResolution.edges.0.source": "losat",
    },
    "visible_text": {
        "#1",
        "Streptomyces lividus CBS 844.73",
        "#2",
        "Streptomyces fradiae ATCC 10745",
        "Run LOSAT",
    },
}


FIRST_RAW_RESULT_CAPTURE = {
    "src": "./media/BGC0000708-BGC0000713/manual-03-04-first-raw-result.webp",
    "selector": "[data-linear-raw-result='record-1->record-2']",
    "viewport": {"width": 1600, "height": 1100},
    "state": {
        "mode": "linear",
        "linearComparisonResolution.edges.0.edgeKey": "record-1->record-2",
        "linearComparisonResolution.edges.0.source": "losat",
    },
    "visible_text": {"Raw LOSAT filename", "Raw result", "Save Raw LOSAT TSV"},
}


def _load_tutorial(sample: dict[str, object]) -> dict[str, object]:
    tutorial_path = resolve_gallery_reference(str(sample["tutorial"]))
    assert tutorial_path is not None
    return json.loads(tutorial_path.read_text(encoding="utf-8"))


def test_linear_comparison_panels_capture_the_command_and_current_status() -> None:
    for example_id in LINEAR_COMPARISON_PANEL_EXAMPLES:
        sample = load_ready_examples(example_id)[0]
        tutorial = _load_tutorial(sample)
        panel_operations = [
            operation
            for _, operation in iter_operation_contexts(sample, tutorial)
            if operation.get("capture", {}).get("selector")
            == "[data-linear-comparison-card]"
            and "Current: Run LOSAT"
            in operation.get("capture", {}).get("visibleText", [])
        ]

        assert len(panel_operations) == 1
        capture = panel_operations[0]["capture"]
        assert panel_operations[0]["dataDependent"] is True
        assert capture["session"] == sample["session"]
        assert {"Comparison", "Run LOSAT", "Current: Run LOSAT", "Settings"} <= set(
            capture["visibleText"]
        )
        assert any(text.startswith("Selected pairs") for text in capture["visibleText"])
        assert "setLinearComparisonGlobalAction('losat')" in " ".join(
            action.get("script", "") for action in capture.get("actions", [])
        )


@pytest.mark.parametrize("example_id", LOSATP_MODE_CAPTURES)
def test_losatp_mode_captures_show_both_levels(
    example_id: str,
) -> None:
    expected = LOSATP_MODE_CAPTURES[example_id]
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    matches = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src") == expected["src"]
    ]

    assert len(matches) == 1
    operation = matches[0]
    capture = operation["capture"]
    scripts = " ".join(
        action.get("script", "") for action in capture.get("actions", [])
    )
    assert capture["crop"] == "openSelect"
    assert capture["cropPadding"]["top"] >= 60
    assert capture["openSelect"] == {
        "selectedLabel": expected["selected"],
        "options": ["Similarity groups", "Collinear blocks", "Pairwise matches"],
        "label": "LOSATP mode",
        "hideFollowingSiblings": True,
    }
    assert capture["clipSelectors"] == ["[data-linear-comparison-losat-mode]"]
    assert capture["assertAppState"]["losatProgram"] == "blastp"
    assert capture["assertAppState"]["losat.blastp.mode"] == expected["state"]
    assert "setLinearComparisonLosatMode('blastp')" in scripts
    assert f"setLinearComparisonLosatpMode('{expected['state']}')" in scripts
    assert capture["visibleControls"] == [
        {
            "selector": "[data-linear-comparison-losat-mode-option='blastn']",
            "pressed": False,
        },
        {
            "selector": "[data-linear-comparison-losat-mode-option='blastp']",
            "pressed": True,
        },
        {
            "selector": "[data-linear-comparison-losat-mode-option='tblastx']",
            "pressed": False,
        },
        {"label": "LOSATP mode", "value": expected["selected"]},
    ]
    assert capture["visibleText"] == ["LOSAT Mode"]
    assert "LOSAT Mode" in operation["media"]["alt"]
    assert "LOSATP mode" in operation["media"]["alt"]


def test_linear_comparison_tutorials_do_not_use_the_retired_search_method_api() -> None:
    for example_id in LINEAR_COMPARISON_PANEL_EXAMPLES:
        sample = load_ready_examples(example_id)[0]
        tutorial = _load_tutorial(sample)
        serialized = json.dumps(tutorial)
        assert "setLinearComparisonSearchMethod" not in serialized
        assert '"label": "Search method"' not in serialized


def test_first_linear_comparison_boundary_has_an_exact_capture_contract() -> None:
    sample = load_ready_examples("BGC0000708-BGC0000713")[0]
    tutorial = _load_tutorial(sample)
    matches = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src")
        == FIRST_COMPARISON_BOUNDARY_CAPTURE["src"]
    ]

    assert len(matches) == 1
    operation = matches[0]
    capture = operation["capture"]
    expected = FIRST_COMPARISON_BOUNDARY_CAPTURE
    assert operation["dataDependent"] is True
    assert capture["source"] == "webapp"
    assert capture["session"] == sample["session"]
    assert capture["selector"] == expected["selector"]
    assert capture["waitForSelector"] == expected["selector"]
    assert capture["viewport"] == expected["viewport"]
    assert capture["cropPadding"] == {
        "top": 0,
        "right": 24,
        "bottom": 12,
        "left": 12,
    }
    assert expected["state"].items() <= capture["assertAppState"].items()
    assert expected["visible_text"] <= set(capture["visibleText"])
    assert capture["openDetails"] == ["Selected pairs"]
    assert capture["visibleControls"] == [
        {
            "selector": (
                "[data-edge-key='record-1->record-2'] "
                "input[type='radio']:checked"
            ),
            "checked": True,
        },
    ]


def test_first_raw_result_has_an_exact_advanced_capture_contract() -> None:
    sample = load_ready_examples("BGC0000708-BGC0000713")[0]
    tutorial = _load_tutorial(sample)
    expected = FIRST_RAW_RESULT_CAPTURE
    matches = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src") == expected["src"]
    ]

    assert len(matches) == 1
    operation = matches[0]
    capture = operation["capture"]
    assert operation["dataDependent"] is True
    assert capture["source"] == "webapp"
    assert capture["session"] == sample["session"]
    assert capture["selector"] == expected["selector"]
    assert capture["waitForSelector"] == expected["selector"]
    assert capture["viewport"] == expected["viewport"]
    assert capture["cropPadding"] == 14
    assert capture["openDetails"] == ["Advanced comparison and layout"]
    assert expected["state"].items() <= capture["assertAppState"].items()
    assert expected["visible_text"] <= set(capture["visibleText"])
    assert capture["visibleControls"] == [
        {
            "selector": (
                "[data-linear-raw-result='record-1->record-2'] "
                "input[aria-label='Raw LOSAT filename for #1 to #2']"
            ),
            "value": "",
        }
    ]


def test_vibrio_capture_keeps_publication_comparisons_active() -> None:
    sample = load_ready_examples("vibrio-harveyi-group-collinear")[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    with gzip.open(session_path, "rt", encoding="utf-8") as handle:
        session = json.load(handle)
    assert session["config"]["linearComparisonPlan"] == {
        "mode": "adjacent",
        "defaultSource": "losat",
        "edges": [],
    }

    tutorial = _load_tutorial(sample)
    comparison_captures = [
        operation["capture"]
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src", "").startswith(
            "./media/vibrio-harveyi-group-collinear/manual-04-"
        )
    ]
    assert comparison_captures
    for capture in comparison_captures:
        scripts = " ".join(
            action.get("script", "") for action in capture.get("actions", [])
        )
        assert "setLinearComparisonGlobalAction('losat')" in scripts
        assert capture["assertAppState"]["linearComparisonPlan.mode"] == "adjacent"


@pytest.mark.parametrize("example_id", LAYOUT_RESULT_CAPTURES)
def test_layout_sensitive_gallery_results_have_exact_capture_contracts(
    example_id: str,
) -> None:
    expected = LAYOUT_RESULT_CAPTURES[example_id]
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    matches = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src") == expected["src"]
    ]

    assert len(matches) == 1
    operation = matches[0]
    capture = operation["capture"]
    assert operation["dataDependent"] is True
    assert capture["source"] == "webapp"
    assert capture["session"] == sample["session"]
    assert capture["selector"] == expected["selector"]
    assert capture["viewport"] == expected["viewport"]
    assert capture["cropPadding"] == 12
    assert expected["state"].items() <= capture["assertAppState"].items()
    assert expected["visible_text"] <= set(capture["visibleText"])


@pytest.mark.parametrize("example_id", BASIC_EXAMPLES)
def test_basic_gallery_operations_have_owned_semantic_capture_contracts(
    example_id: str,
) -> None:
    expected = BASIC_EXAMPLES[example_id]
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    operations = [
        operation
        for _, operation in iter_operation_contexts(sample, tutorial)
    ]
    data_dependent = [
        operation for operation in operations if operation.get("dataDependent") is True
    ]

    assert len(operations) == expected["operation_count"]
    assert len(data_dependent) == expected["data_dependent_count"]

    for operation in data_dependent:
        media = operation["media"]
        capture = operation["capture"]
        assert media["src"].startswith(f"./media/{example_id}/")
        assert capture["source"] == "webapp"
        assert capture["session"] == sample["session"]
        assert capture["assertAppState"]
        assert capture.get("visibleControls") or capture.get("visibleText")

    visible_identity = {
        text
        for operation in data_dependent
        for text in operation["capture"].get("visibleText", [])
    }
    assert expected["filename"] in visible_identity


@pytest.mark.parametrize("example_id", BASIC_EXAMPLES)
def test_basic_gallery_capture_session_contains_the_named_record(
    example_id: str,
) -> None:
    expected = BASIC_EXAMPLES[example_id]
    sample = load_ready_examples(example_id)[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    session = json.loads(session_path.read_text(encoding="utf-8"))
    resource = session["resources"]["record-1-genbank"]

    assert resource["encoding"] == "base64"
    decoded = base64.b64decode(resource["data"])
    assert expected["resource_marker"] in decoded


def test_ready_gallery_media_passes_capture_contract_validation() -> None:
    result = validate_tutorial_media(
        load_ready_examples(),
        max_width=DEFAULT_MAX_IMAGE_WIDTH,
        max_height=DEFAULT_MAX_IMAGE_HEIGHT,
        max_file_size_kb=DEFAULT_MAX_FILE_SIZE_KB,
    )

    assert result.errors == []


def test_data_dependent_capture_rejects_cross_example_media() -> None:
    context = TutorialContext(
        example_id="current-example",
        tutorial_path=Path("tutorial.json"),
        section="manualSteps",
        step_index=1,
        step_title="Upload",
        operation_index=1,
    )
    sample = {"session": "./sessions/current-example.gbdraw-session.json"}
    operation = {
        "dataDependent": True,
        "media": {"src": "./media/other-example/upload.webp"},
        "capture": {
            "source": "webapp",
            "session": "./sessions/current-example.gbdraw-session.json",
            "assertAppState": {"files.c_gb.name": "current.gbk"},
            "visibleText": ["current.gbk"],
        },
    }
    result = ValidationResult(errors=[], warnings=[])

    add_capture_contract_validation(result, context, sample, operation)

    assert any(
        "data-dependent media must be owned by current-example" in error
        for error in result.errors
    )


def test_cross_example_media_is_explicitly_generic() -> None:
    for sample in load_ready_examples():
        tutorial = _load_tutorial(sample)
        for _, operation in iter_operation_contexts(sample, tutorial):
            media = operation.get("media")
            if not isinstance(media, dict):
                continue
            source = str(media.get("src", ""))
            if not source.startswith("./media/"):
                continue
            media_owner = source.split("/", 3)[2]
            if media_owner != sample["id"]:
                assert operation.get("genericMedia") is True


def test_web_app_readiness_waits_for_palette_without_starting_worker() -> None:
    calls: list[tuple[str, int]] = []

    class FakePage:
        def wait_for_function(self, script: str, *, timeout: int) -> None:
            calls.append((script, timeout))

    wait_for_app_shell(FakePage())

    assert len(calls) == 1
    script, timeout = calls[0]
    assert "Object.keys(app.paletteDefinitions || {}).length > 0" in script
    assert "diagramGenerationWorkerReady" not in script
    assert "pyodideReady" not in script
    assert timeout == 120000
