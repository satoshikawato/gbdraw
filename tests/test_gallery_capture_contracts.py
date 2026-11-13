from __future__ import annotations

import base64
import io
import json
from pathlib import Path

import pytest
from PIL import Image

import tools.capture_gallery_tutorial_screenshots as capture_tool
from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import read_session_payload
from tools.capture_gallery_tutorial_screenshots import (
    DEFAULT_MAX_FILE_SIZE_KB,
    DEFAULT_MAX_IMAGE_HEIGHT,
    DEFAULT_MAX_IMAGE_WIDTH,
    TutorialContext,
    ValidationResult,
    add_capture_contract_validation,
    add_step_replay_validation,
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
        "src": "./media/tobacco-chloroplast/manual-08-02-chloroplast-preview.webp",
        "selector": "svg[data-capture-chloroplast-preview='true']",
        "viewport": {"width": 1400, "height": 1400},
        "state": {
            "mode": "circular",
            "form.legend": "upper_left",
            "annotationSets.0.id": "plastome_regions",
            "adv.circular_track_slots_enabled": True,
        },
        "visible_text": {"Nicotiana tabacum", "NC_001879.2", "LSC", "IRb", "SSC", "IRa", "GC content"},
    },
    "Vnig_TUMSAT-TG-2018": {
        "src": "./media/Vnig_TUMSAT-TG-2018/manual-06-02-multirecord-preview.webp",
        "selector": "svg[data-capture-vnig-preview='true']",
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
                {"selector": "record-1:1", "row": 1},
                {"selector": "record-1:2", "row": 1},
                {"selector": "record-1:3", "row": 2},
                {"selector": "record-1:4", "row": 2},
                {"selector": "record-1:5", "row": 2},
                {"selector": "record-1:6", "row": 2},
            ],
            "form.plot_title": "<i>Vibrio nigripulchritudo</i> TUMSAT-TG-2018, complete genome",
            "adv.plot_title_position": "bottom",
        },
        "visible_text": {"Vibrio nigripulchritudo", "AP024087.1", "AP024088.1", "pVNTG1", "pVNTG4", "GC content", "GC skew (+)"},
    },
    "vibrio-harveyi-group-collinear": {
        "src": "./media/vibrio-harveyi-group-collinear/manual-10-01-collinear-overview.webp",
        "selector": "svg[data-capture-vibrio-preview='true']",
        "viewport": {"width": 1800, "height": 1100},
        "state": {
            "mode": "linear",
            "linearSeqs.length": 4,
            "linearRecordLayoutEnabled": True,
            "linearRecordGap": 48,
            "linearRecordRows.length": 4,
            "form.legend": "bottom",
            "form.plot_title": "",
            "form.linear_track_layout": "above",
            "form.scale_style": "ruler",
            "form.show_labels_linear": "none",
            "losatProgram": "blastp",
            "losat.blastp.mode": "collinear",
            "losat.blastp.collinearColorMode": "orientation_identity",
        },
        "visible_text": {"Vibrio parahaemolyticus", "Vibrio alginolyticus", "Collinear", "Inverted"},
    },
    "hepatoplasmataceae_collinear": {
        "src": "./media/hepatoplasmataceae_collinear/manual-07-03-collinear-overview.webp",
        "selector": "svg[data-capture-tutorial-overview='true']",
        "viewport": {"width": 1900, "height": 1000},
        "state": {
            "mode": "linear",
            "linearSeqs.length": 5,
            "form.legend": "right",
            "form.show_gc": True,
            "form.show_skew": True,
            "form.scale_style": "ruler",
            "form.align_center": True,
            "selectedPalette": "ajisai",
            "losat.blastp.mode": "collinear",
            "losat.blastp.collinearColorMode": "orientation_identity",
        },
        "visible_text": {"Tyloplasma litorale", "AP027078.1", "NZ_CP006932.1", "Collinear", "Inverted"},
    },
    "hepatoplasma_scabrum_depth": {
        "src": "./media/hepatoplasma_scabrum_depth/manual-07-02-final-preview.webp",
        "selector": "svg[data-capture-depth-preview='true']",
        "viewport": {"width": 1400, "height": 1400},
        "state": {
            "mode": "circular",
            "form.show_depth": True,
            "adv.depth_max": 80,
            "adv.depth_width_circular": 60,
        },
        "visible_text": {"AP027133.1", "DRR394922 mean depth"},
    },
    "ecoli_o157_sakai_rings": {
        "src": "./media/ecoli_o157_sakai_rings/manual-06-03-final-preview.webp",
        "selector": "svg[data-capture-rings-preview='true']",
        "viewport": {"width": 1400, "height": 1400},
        "state": {
            "mode": "circular",
            "adv.circular_track_slots_enabled": True,
        },
        "visible_text": {"O157:H7 EDL933", "K-12 MG1655", "CFT073 (UPEC)", "Stx2 phage", "LEE island"},
    },
    "HmmtDNA_ATskew": {
        "src": "./media/HmmtDNA_ATskew/manual-09-02-atskew-preview.webp",
        "selector": "svg[data-capture-hmmt-preview='true']",
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
        "src": "./media/BGC0000708-BGC0000713/manual-05-02-select-losatp-orthogroups.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "hepatoplasmataceae_orthogroup": {
        "src": "./media/hepatoplasmataceae_orthogroup/manual-04-01-orthogroups-mode.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "hepatoplasmataceae_collinear": {
        "src": "./media/hepatoplasmataceae_collinear/manual-05-02-collinear-reduction.webp",
        "selected": "Collinear blocks",
        "state": "collinear",
    },
    "majanivirus_orthogroup": {
        "src": "./media/majanivirus_orthogroup/manual-03-02-losatp-orthogroups.webp",
        "selected": "Similarity groups",
        "state": "orthogroup",
    },
    "vibrio-harveyi-group-collinear": {
        "src": "./media/vibrio-harveyi-group-collinear/manual-05-02-search-method-collinear.webp",
        "selected": "Collinear blocks",
        "state": "collinear",
    },
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
            and "setLinearComparisonGlobalAction('losat')" in " ".join(
                action.get("script", "") for action in operation.get("capture", {}).get("actions", [])
            )
        ]

        assert len(panel_operations) == 1
        capture = panel_operations[0]["capture"]
        assert panel_operations[0]["dataDependent"] is True
        assert capture["session"] == sample["session"]
        assert {"Comparison", "Run LOSAT", "Settings"} <= set(
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


def test_bgc_previews_replay_the_reader_steps_through_their_own_generate() -> None:
    """D-01/D-08: each BGC preview is the result of the reader's own steps from the MIBiG files."""

    sample = load_ready_examples("BGC0000708-BGC0000713")[0]
    tutorial = _load_tutorial(sample)
    steps = tutorial["manualSteps"]
    generate_steps = set()
    for number, step in enumerate(steps, start=1):
        types = [action["type"] for action in step.get("replay", [])]
        if "generate" in types:
            generate_steps.add(number)
            assert types[-1] == "generate", number
    # Step 7 regenerates through Apply in the Review alignment dialog, not Generate Diagram.
    assert generate_steps == {4, 5, 6, 9, 10, 11}
    assert any("similarity-alignment-apply-reason" in action.get("selector", "") for action in steps[6]["replay"])

    previews = {
        number: operation["capture"]
        for number, step in enumerate(steps, start=1)
        for operation in step.get("operations", [])
        if operation.get("capture", {}).get("selector") == "[data-capture='bgc-preview']"
    }
    assert set(previews) == generate_steps | {7}
    for number, capture in previews.items():
        assert capture["replayThrough"] == number
        assert capture["session"] == sample["session"]
        assert capture["assertAppState"]["linearSeqs.length"] == 5
    assert previews[5]["assertAppState"]["linearSeqs.4.region_reverse"] is False
    assert previews[6]["assertAppState"]["linearSeqs.4.region_reverse"] is True


def test_bgc_session_holds_the_mibig_files_with_record_reverse_complement() -> None:
    """D-08: the Session stores the original BGC0000713.gbk and reverses it per record."""

    sample = load_ready_examples("BGC0000708-BGC0000713")[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    session = drawing_view(read_session_payload(session_path))
    records = session["renderRequest"]["records"]
    assert [record["presentation"]["reverseComplement"] for record in records] == [
        False,
        False,
        False,
        False,
        True,
    ]
    data_dir = Path("gbdraw/web/tutorial-data/aminoglycoside-bgc-five")
    for record, name in zip(records, sorted(path.name for path in data_dir.glob("*.gbk")), strict=True):
        resource = session["resources"][record["source"]["resourceId"]]
        assert resource["name"].endswith(name), (resource["name"], name)
        assert base64.b64decode(resource["data"]) == (data_dir / name).read_bytes(), name


def test_bgc_finish_step_types_the_arrow_geometry_of_the_published_figure() -> None:
    """The published figure fixes both arrow ratios, so the reader's last step types both.

    Left at Auto, the head length follows the size class of the longest cluster.
    """

    sample = load_ready_examples("BGC0000708-BGC0000713")[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    session = drawing_view(read_session_payload(session_path))
    published = session["renderRequest"]["diagramOptions"]["configOverrides"]
    step = _load_tutorial(sample)["manualSteps"][10]
    typed = {action["selector"]: action["value"] for action in step["replay"] if action["type"] == "fill"}
    finished = step["operations"][-1]
    for name, label in (("head_length_ratio", "Head Length Ratio"), ("shaft_width_ratio", "Shaft Width Ratio")):
        value = published[f"objects.features.arrow_geometry.{name}"]
        row = [f"Features › Arrow Geometry › {label}", str(value)]
        assert typed[f'input[aria-label="Arrow {label.lower()}"]'] == value
        assert row in step["table"]["rows"]
        assert row in finished["table"]["rows"]
        assert finished["capture"]["assertAppState"][f"adv.arrow_{name}"] == value


def test_vibrio_capture_keeps_publication_comparisons_active() -> None:
    sample = load_ready_examples("vibrio-harveyi-group-collinear")[0]
    session_path = resolve_gallery_reference(str(sample["session"]))
    assert session_path is not None
    session = drawing_view(read_session_payload(session_path))
    assert session["modes"]["linear"]["config"]["linearComparisonPlan"] == {
        "mode": "adjacent",
        "defaultSource": "losat",
        "edges": [],
    }

    tutorial = _load_tutorial(sample)
    comparison_captures = [
        operation["capture"]
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src", "").startswith(
            "./media/vibrio-harveyi-group-collinear/manual-05-"
        )
        and "replayThrough" not in operation["capture"]
    ]
    assert len(comparison_captures) == 4
    for capture in comparison_captures:
        scripts = " ".join(
            action.get("script", "") for action in capture.get("actions", [])
        )
        assert "setLinearComparisonGlobalAction('losat')" in scripts
        assert capture["assertAppState"]["linearComparisonPlan.mode"] == "adjacent"


def test_vibrio_previews_replay_each_generate_step() -> None:
    """Every Generate step shows its own replayed preview; LOSAT is set up once (D-01, D-02, D-07)."""

    sample = load_ready_examples("vibrio-harveyi-group-collinear")[0]
    tutorial = _load_tutorial(sample)
    steps = tutorial["manualSteps"]
    generate_steps = [
        index
        for index, step in enumerate(steps, start=1)
        if any(action["type"] == "generate" for action in step.get("replay", []))
    ]
    assert generate_steps == [4, 5, 8, 10]

    losat_steps = [
        index
        for index, step in enumerate(steps, start=1)
        if any("Run LOSAT" in action.get("selector", "") for action in step.get("replay", []))
    ]
    assert losat_steps == [5]

    # File defaults are filled one File at a time; no step adds a slot during replay.
    assert [action["selector"] for action in steps[1]["replay"]] == [
        "input[aria-label='Default definition for file 1']",
        "input[aria-label='Default subtitle for file 1']",
    ]
    assert [action["selector"] for action in steps[2]["replay"]] == [
        "input[aria-label='Default definition for file 2']",
        "input[aria-label='Default subtitle for file 2']",
    ]

    for index in generate_steps:
        previews = [
            operation["capture"]
            for operation in steps[index - 1]["operations"]
            if operation["capture"].get("selector") == "svg[data-capture-vibrio-preview='true']"
        ]
        assert len(previews) == 1, index
        capture = previews[0]
        assert capture["replayThrough"] == index
        assert capture["session"] == sample["session"]
        scripts = " ".join(action.get("script", "") for action in capture["actions"])
        assert "(D-07)" in scripts


PROGRESSIVE_TUTORIALS = {
    # Generate steps, the preview selector, and the replay selectors the finished figure needs
    # beyond the Web defaults (OV-1695, OV-1699: Multi-Record Canvas off).
    "hepatoplasmataceae_collinear": {
        "generate_steps": [4, 5, 6, 7],
        "selector": "svg[data-capture-tutorial-overview='true']",
        # LOSAT settings set before the first run, so one LOSAT run serves every later Generate.
        "required_selectors": {
            "button[aria-label='Run LOSAT for all adjacent pairs']",
            "input[aria-label='LOSATP Max target seqs']",
            "input[aria-label='Build similarity groups with within-record searches']",
            "input[aria-label='Lock Definition Column']",
        },
    },
    "Vnig_TUMSAT-TG-2018": {
        "generate_steps": [3, 4, 5, 6],
        "selector": "svg[data-capture-vnig-preview='true']",
        "required_selectors": {"select[aria-label^='Row for #3 (']"},
    },
    "HmmtDNA_ATskew": {
        "generate_steps": [3, 4, 5, 6, 7, 9],
        "selector": "svg[data-capture-hmmt-preview='true']",
        "required_selectors": {"input[data-circular-canvas-setting]", "input[aria-label='Separate Strands']"},
    },
    "hepatoplasma_scabrum_depth": {
        "generate_steps": [3, 4, 5, 6, 7],
        "selector": "svg[data-capture-depth-preview='true']",
        # The depth file stays through Reset with Show Depth off; the depth step chooses
        # it again, as the reader does, and turns Show Depth on.
        "required_selectors": {
            "role=button[name='Choose Depth TSV'] >> visible=true",
            "role=checkbox[name='Show Depth'] >> visible=true",
            "input[aria-label='Depth Window'] >> visible=true",
            "input[aria-label='Depth large tick']",
        },
    },
    "ecoli_o157_sakai_rings": {
        "generate_steps": [3, 4, 5, 6],
        "selector": "svg[data-capture-rings-preview='true']",
        # Reset keeps the comparison files with Show Comparison Rings off and the ring names back to
        # the accessions; it removes the imported region sets and turns the custom stack off.
        "required_selectors": {
            "role=checkbox[name='Show Comparison Rings'] >> visible=true",
            "input[aria-label='Comparison ring label 1']",
            "button:has-text('Import TSV')",
            "[aria-label='Circular track slot annotations_1 Extend to']",
        },
    },
    "tobacco-chloroplast": {
        "generate_steps": [3, 4, 5, 6, 7, 8],
        "selector": "svg[data-capture-chloroplast-preview='true']",
        "required_selectors": {
            "input[data-circular-canvas-setting]",
            "button[aria-label='Remove repeat_region']",
            "input[aria-label='Block stroke color']",
            "input[aria-label='Circular Label Spacing']",
        },
    },
}


@pytest.mark.parametrize("example_id", PROGRESSIVE_TUTORIALS)
def test_progressive_tutorials_capture_each_generate_by_replay(example_id: str) -> None:
    """Each Generate step shows its own preview, replayed from the Web defaults (TV-D18)."""

    expected = PROGRESSIVE_TUTORIALS[example_id]
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    steps = tutorial["manualSteps"]
    generate_steps = []
    for number, step in enumerate(steps, start=1):
        types = [action["type"] for action in step.get("replay", [])]
        if "generate" in types:
            generate_steps.append(number)
            assert types.count("generate") == 1 and types[-1] == "generate", number
    assert generate_steps == expected["generate_steps"]

    for number in generate_steps:
        previews = [
            operation
            for operation in steps[number - 1]["operations"]
            if operation.get("capture", {}).get("selector") == expected["selector"]
        ]
        assert len(previews) == 1, number
        capture = previews[0]["capture"]
        assert previews[0]["dataDependent"] is True
        assert capture["replayThrough"] == number
        assert capture["session"] == sample["session"]
    final = steps[generate_steps[-1] - 1]["operations"][-1]
    assert final["media"]["src"] == LAYOUT_RESULT_CAPTURES[example_id]["src"]

    replayed = [action for step in steps for action in step.get("replay", [])]
    assert expected["required_selectors"] <= {action.get("selector") for action in replayed}
    downloads = {entry["targetName"] for entry in tutorial["downloads"]}
    for action in replayed:
        if action["type"] == "upload":
            assert action["path"].startswith("gallery/files/")
            assert Path(action["path"]).name in downloads


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
    session = drawing_view(read_session_payload(session_path))
    resource = session["resources"]["record-1-genbank"]

    assert resource["encoding"] == "base64"
    decoded = base64.b64decode(resource["data"])
    assert expected["resource_marker"] in decoded


def test_vnig_upload_capture_names_the_session_file_it_shows() -> None:
    """The tool retakes the upload crop only from capture metadata (OV-628)."""

    filename = "GCF_015097735.1_ASM1509773v1_genomic.gbff"
    sample = load_ready_examples("Vnig_TUMSAT-TG-2018")[0]
    operations = [
        operation
        for _, operation in iter_operation_contexts(sample, _load_tutorial(sample))
        if operation.get("media", {}).get("src", "").endswith("/manual-02-01-gbff-upload.webp")
    ]

    assert len(operations) == 1
    operation = operations[0]
    assert operation.get("dataDependent") is True
    capture = operation.get("capture", {})
    assert capture.get("source") == "webapp"
    assert capture.get("session") == sample["session"]
    # The file of the first Circular File card (OV-3328: `files.c_gb` is gone), and the six
    # records the step says the app finds, read as the reader's upload reads them (OV-3139).
    assert capture.get("assertAppState") == {
        "mode": "circular",
        "circularFileCards.cards.value.0.source.gb.name": filename,
    }
    assert capture.get("visibleText") == [filename, "6 source records inspected"]


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
    sample = {"session": "./sessions/current-example.gbdraw"}
    operation = {
        "dataDependent": True,
        "media": {"src": "./media/other-example/upload.webp"},
        "capture": {
            "source": "webapp",
            "session": "./sessions/current-example.gbdraw",
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


def _manual_context(step_index: int) -> TutorialContext:
    return TutorialContext(
        example_id="current-example",
        tutorial_path=Path("tutorial.json"),
        section="manualSteps",
        step_index=step_index,
        step_title="Generate",
        operation_index=1,
    )


@pytest.mark.parametrize("replay_through", [0, 3, "2", True])
def test_replay_capture_names_a_step_up_to_its_own(replay_through: object) -> None:
    result = ValidationResult(errors=[], warnings=[])
    operation = {"capture": {"source": "webapp", "replayThrough": replay_through}}

    add_capture_contract_validation(result, _manual_context(2), {}, operation)

    assert result.errors == [
        "current-example / manualSteps / step 2 / Generate / operation 1: "
        "capture.replayThrough must name a manual step from 1 to this step."
    ]


def test_step_replay_accepts_reader_actions_only() -> None:
    result = ValidationResult(errors=[], warnings=[])
    step = {
        "replay": [
            {"type": "fill"},
            {"type": "generate"},
            {"type": "resetSettings"},
            {"type": "upload", "selector": "role=button[name='Import TSV']", "path": "gallery/files/qualifier_priority.tsv"},
            {"type": "upload", "selector": "role=button[name='Import TSV']"},
            {"type": "upload", "selector": "role=button[name='Import TSV']", "path": "../../pyproject.toml"},
            {"type": "upload", "selector": "role=button[name='Import TSV']", "path": "gallery/files/missing.tsv"},
            {"type": "upload", "path": "gallery/files/qualifier_priority.tsv"},
        ]
    }

    add_step_replay_validation(result, _manual_context(1), step)

    assert [error.split(": ", 1)[1].split(" ")[2] for error in result.errors] == ["3", "5", "6", "7", "8"]
    assert "unsupported type" in result.errors[0]
    assert all("upload" in error for error in result.errors[1:])


def test_upload_action_chooses_the_file_through_the_named_button() -> None:
    """A reader's table upload: click the button that opens the file chooser, then choose the file."""

    events: list[object] = []

    class FakeChooser:
        def set_files(self, path: str) -> None:
            events.append(("set_files", path))

    class FakeChooserInfo:
        value = FakeChooser()

        def __enter__(self) -> FakeChooserInfo:
            events.append("expect_file_chooser")
            return self

        def __exit__(self, *exc: object) -> None:
            events.append("file_chooser_opened")

    class FakeLocator:
        def __init__(self, selector: str) -> None:
            self.selector = selector

        @property
        def first(self) -> FakeLocator:
            return self

        def click(self) -> None:
            events.append(("click", self.selector))

    class FakePage:
        def expect_file_chooser(self) -> FakeChooserInfo:
            return FakeChooserInfo()

        def locator(self, selector: str) -> FakeLocator:
            return FakeLocator(selector)

        def evaluate(self, script: str) -> None:
            events.append("wait for transitions")

    selector = "role=button[name='Choose Priority File (TSV)']"
    capture_tool.apply_capture_action(
        FakePage(), {"type": "upload", "selector": selector, "path": "gallery/files/qualifier_priority.tsv"}
    )

    expected = capture_tool.WEB_ROOT.resolve() / "gallery/files/qualifier_priority.tsv"
    assert events == [
        "expect_file_chooser",
        ("click", selector),
        "file_chooser_opened",
        ("set_files", str(expected)),
        "wait for transitions",
    ]


def test_replay_resets_then_performs_steps_through_the_named_one(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    tutorial = tmp_path / "tutorial.json"
    tutorial.write_text(
        json.dumps(
            {
                "manualSteps": [
                    {"replay": [{"type": "fill", "value": "one"}]},
                    {"title": "no replay"},
                    {"replay": [{"type": "generate"}]},
                    {"replay": [{"type": "fill", "value": "four"}]},
                ]
            }
        )
    )
    performed: list[dict[str, object]] = []
    monkeypatch.setattr(capture_tool, "resolve_gallery_reference", lambda reference: tutorial)
    monkeypatch.setattr(capture_tool, "apply_capture_action", lambda page, action: performed.append(action))

    capture_tool.replay_tutorial_steps(object(), {"tutorial": "./tutorials/x.json"}, 3)

    assert performed == [
        {"type": "resetSettings"},
        {"type": "fill", "value": "one"},
        {"type": "generate"},
    ]


def test_reset_settings_waits_for_the_diagram_update_to_finish() -> None:
    """Reset Settings redraws the shown Result; the mode buttons and Undo stay disabled until it ends."""

    events: list[str] = []

    class FakeLocator:
        def __init__(self, name: str) -> None:
            self.name = name

        def click(self) -> None:
            events.append(f"click {self.name}")

        def get_by_role(self, role: str, name: str, exact: bool = False) -> FakeLocator:
            return FakeLocator(name)

        def wait_for(self, state: str) -> None:
            events.append(f"wait_for {state}")

    class FakePage:
        def get_by_role(self, role: str, name: str, exact: bool = False) -> FakeLocator:
            return FakeLocator(name)

        def wait_for_timeout(self, ms: int) -> None:
            events.append(f"timeout {ms}")

        def wait_for_function(self, script: str, timeout: int) -> None:
            events.append(script)

        def evaluate(self, script: str) -> None:
            events.append(script)

    capture_tool.apply_capture_action(FakePage(), {"type": "resetSettings"})

    assert events[:3] == ["click Reset Settings", "click Reset all drawings", "wait_for hidden"]
    waits = [event for event in events[3:] if "labelReflowProcessing" in event]
    assert len(waits) == 1
    assert "processing" in waits[0].replace("labelReflowProcessing", "")
    # After the redraw, History can stay busy for many seconds (OV-1700); the replay's next
    # step needs the drawing tabs, so Reset also waits until they are enabled.
    switch_waits = [event for event in events[3:] if ".drawing-bar [role=" in event]
    assert len(switch_waits) == 1 and "disabled" in switch_waits[0]


def test_vibrio_rotation_steps_follow_the_open_popup_and_close_it() -> None:
    """With the popup open, > retargets it to the next match; the reader closes it before the next step."""

    sample = load_ready_examples("vibrio-harveyi-group-collinear")[0]
    steps = _load_tutorial(sample)["manualSteps"]
    for number in (6, 7):
        replay = steps[number - 1]["replay"]
        selectors = [action.get("selector", "") for action in replay]
        next_index = selectors.index("button[aria-label='Next match']")
        assert "button[aria-label='Open active feature']" in selectors[:next_index], number
        assert "button[aria-label='Open active feature']" not in selectors[next_index:], number
        close_index = selectors.index("button[aria-label='Close feature popup']")
        assert close_index > next_index, number
        assert replay[close_index + 1] == {
            "type": "waitForSelector",
            "selector": ".feature-popup",
            "state": "hidden",
        }, number
        assert "close the popup" in steps[number - 1]["body"], number
    for operation in steps[6]["operations"]:
        scripts = " ".join(action.get("script", "") for action in operation["capture"].get("actions", []))
        assert "Close feature popup" not in scripts


def test_capture_lets_css_transitions_end_before_the_crop(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A caret turns over 200 ms after its card opens; a crop taken sooner shows it mid-turn (OV-627)."""

    events: list[str] = []

    class Mouse:
        def move(self, x: float, y: float) -> None:
            pass

    class FakePage:
        mouse = Mouse()

        def evaluate(self, script: str, *args: object) -> list[object]:
            events.append(script)
            return []

        def screenshot(self, **kwargs: object) -> bytes:
            events.append("screenshot")
            buffer = io.BytesIO()
            Image.new("RGB", (2, 2)).save(buffer, "PNG")
            return buffer.getvalue()

    monkeypatch.setattr(capture_tool, "prepare_capture_page", lambda *args: None)
    monkeypatch.setattr(capture_tool, "page_scroll_dimensions", lambda page: {"width": 2.0, "height": 2.0})
    monkeypatch.setattr(capture_tool, "output_path_for_operation", lambda operation: tmp_path / "crop.webp")
    operation = {"capture": {"boundingBox": {"x": 0, "y": 0, "width": 2, "height": 2}}}

    capture_tool.capture_operation(FakePage(), {}, operation, "http://127.0.0.1", 90)

    waits = [
        index
        for index, script in enumerate(events)
        if "CSSTransition" in script and ".finished" in script
    ]
    assert waits and waits[-1] < events.index("screenshot")


def test_vnig_palette_select_reads_the_app_palette_list() -> None:
    tutorial = json.loads(
        (Path(__file__).resolve().parents[1] / "gbdraw/web/gallery/tutorials/Vnig_TUMSAT-TG-2018.json").read_text()
    )
    specs = [
        operation["capture"]["openSelect"]
        for step in tutorial["manualSteps"]
        for operation in step.get("operations", [])
        if isinstance(operation, dict)
        and isinstance(operation.get("capture", {}).get("openSelect"), dict)
        and operation["capture"]["openSelect"].get("selectedLabel") == "orchid"
    ]

    assert len(specs) == 1
    assert "options" not in specs[0]


@pytest.mark.parametrize(
    "example_id",
    [
        "BGC0000708-BGC0000713",
        "hepatoplasmataceae_collinear",
        "hepatoplasmataceae_orthogroup",
        "majanivirus_orthogroup",
    ],
)
def test_linear_drawing_bar_crop_is_taken_from_a_new_page_after_the_plus_linear_click(
    example_id: str,
) -> None:
    sample = load_ready_examples(example_id)[0]
    tutorial = _load_tutorial(sample)
    captures = [
        operation["capture"]
        for _, operation in iter_operation_contexts(sample, tutorial)
        if operation.get("media", {}).get("src", "").endswith("manual-01-01-linear-mode.webp")
    ]
    assert len(captures) == 1
    capture = captures[0]
    assert capture["freshPage"] is True
    assert "session" not in capture and "replayThrough" not in capture
    assert capture["actions"][0]["selector"].endswith('[data-drawing-offer="linear"]')


def test_capture_clip_refuses_a_crop_that_cuts_a_clipped_element() -> None:
    """A popup placed past the page edge was written cut at the edge (OV-2280); the tool refuses it."""

    padding = {"top": 18, "right": 180, "bottom": 2, "left": 18}
    page = {"width": 1280.0, "height": 900.0}
    inside = capture_tool.padded_clip_from_boxes([{"x": 50, "y": 255, "width": 566, "height": 62}], padding, page)
    assert inside == {"x": 32.0, "y": 237.0, "width": 764.0, "height": 82.0}
    # Padding past the page edge is trimmed; only the element itself must fit.
    edge = capture_tool.padded_clip_from_boxes([{"x": 1000, "y": 10, "width": 280, "height": 20}], padding, page)
    assert edge["x"] + edge["width"] == 1280.0
    with pytest.raises(RuntimeError, match="cut"):
        capture_tool.padded_clip_from_boxes([{"x": 1010, "y": 255, "width": 566, "height": 62}], padding, page)


def test_popup_captures_find_popup_parts_by_their_hooks() -> None:
    """A capture finds the part of a popup it crops by its hook, never by its place among siblings."""

    positional = ("children[", "firstElementChild", "lastElementChild")
    tutorials = sorted((Path(__file__).resolve().parents[1] / "gbdraw/web/gallery/tutorials").glob("*.json"))
    found = []

    def scan(node: object, where: str) -> None:
        if isinstance(node, dict):
            script = node.get("script")
            if isinstance(script, str):
                found.extend(f"{where}: {name}" for name in positional if name in script)
            for key, value in node.items():
                scan(value, f"{where}.{key}")
        elif isinstance(node, list):
            for index, value in enumerate(node):
                scan(value, f"{where}[{index}]")

    for path in tutorials:
        scan(json.loads(path.read_text(encoding="utf-8")), path.name)
    assert tutorials and not found, found


def test_each_capture_action_lets_its_transitions_end_before_the_next() -> None:
    """A zoom or pan moves the preview in a 0.2 s transition. The popup script measured the
    ribbon before it ended, at its old place off the page, and put the popup past the edge (OV-2280)."""

    events: list[str] = []

    class FakePage:
        def evaluate(self, script: str, *args: object) -> None:
            events.append(script)

    capture_tool.apply_capture_action(FakePage(), {"type": "evaluate", "script": "app.zoom = 1.15"})

    assert events[0] == "app.zoom = 1.15"
    assert len(events) == 2 and "CSSTransition" in events[1] and ".finished" in events[1]


class _PreparePage:
    """A page that records what prepare_capture_page does to it, in order."""

    def __init__(self, events: list[str]) -> None:
        self.events = events

    def set_viewport_size(self, size: dict[str, int]) -> None:
        pass

    def goto(self, url: str, wait_until: str) -> None:
        self.events.append("open page")

    def evaluate(self, script: str, *args: object) -> None:
        if "data-linear-source-fold" in script:
            self.events.append("open Linear File cards")

    def add_style_tag(self, content: str) -> None:
        if "[data-session-file-picker]" in content and "display: none" in content:
            self.events.append("hide Files in this Session")


@pytest.mark.parametrize("fresh_page", [False, True])
def test_session_capture_reads_the_records_as_the_readers_upload_did(
    fresh_page: bool, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-3139: a loaded Session leaves its records unread; the reader who uploaded the files sees them read."""

    events: list[str] = []
    monkeypatch.setattr(capture_tool, "wait_for_app_shell", lambda page: None)
    monkeypatch.setattr(capture_tool, "load_web_app_session", lambda page, path: events.append("load session"))
    monkeypatch.setattr(capture_tool, "dismiss_session_notice", lambda page: events.append("dismiss Session notice"))
    monkeypatch.setattr(capture_tool, "inspect_shown_records", lambda page: events.append("read records"))
    monkeypatch.setattr(capture_tool, "replay_tutorial_steps", lambda page, sample, through: events.append("replay"))
    monkeypatch.setattr(capture_tool, "apply_capture_action", lambda page, action: events.append("action"))
    sample = load_ready_examples("HmmtDNA_ATskew")[0]
    capture: dict[str, object] = {"source": "webapp", "actions": [{"type": "wait"}]}
    if fresh_page:
        capture["freshPage"] = True
    else:
        capture.update(session=sample["session"], replayThrough=2)

    capture_tool.prepare_capture_page(_PreparePage(events), "http://gallery", sample, capture)

    if fresh_page:
        assert events == ["open page", "action"]
    else:
        assert events == [
            "open page",
            "load session",
            "dismiss Session notice",
            "open Linear File cards",
            "hide Files in this Session",
            "replay",
            "read records",
            "action",
        ]


def test_session_notice_is_closed_by_its_dismiss_button() -> None:
    """OV-3369: the reader who uploaded the files sees no `Session loaded` notice above the Result.

    The notice moved the Preview down under popups that captures place at a fixed height, so their
    crops took a strip of the Preview toolbar.
    """

    events: list[str] = []

    class Locator:
        def __init__(self, selector: str) -> None:
            self.selector = selector

        def locator(self, selector: str) -> "Locator":
            return Locator(f"{self.selector} {selector}")

        def click(self) -> None:
            events.append(f"click {self.selector}")

        def wait_for(self, state: str) -> None:
            events.append(f"wait {self.selector} {state}")

    class FakePage:
        def locator(self, selector: str) -> Locator:
            return Locator(selector)

    capture_tool.dismiss_session_notice(FakePage())

    assert events == ["click [data-app-notice] .app-notice-dismiss", "wait [data-app-notice] hidden"]


def test_webapp_crop_refuses_records_not_inspected() -> None:
    """The guard reads the shown drawing's records in both modes and refuses a crop that shows them unread."""

    class FakePage:
        def __init__(self, read: bool) -> None:
            self.read = read
            self.scripts: list[str] = []

        def evaluate(self, script: str) -> bool:
            self.scripts.append(script)
            return self.read

    unread = FakePage(read=False)
    with pytest.raises(RuntimeError, match="records not inspected"):
        capture_tool.assert_capture_records_read(unread)
    capture_tool.assert_capture_records_read(FakePage(read=True))

    script = unread.scripts[0]
    assert "circularRecordDiscoveryState" in script and "status === 'ready'" in script
    assert "linearRecordSelectorDisabled" in script


def test_no_tutorial_step_names_the_session_file_chooser_that_captures_hide() -> None:
    """A capture stands in for a reader who uploaded the files, who has no Session: it hides
    **Files in this Session** (OV-3330), so no step may tell the reader to use it."""

    result = ValidationResult(errors=[], warnings=[])
    capture_tool.add_step_text_validation(
        result, _manual_context(2), {"title": "Upload", "body": "Pick it in Files in this Session."}
    )
    capture_tool.add_step_text_validation(result, _manual_context(3), {"title": "Upload", "body": "Add the file."})

    assert len(result.errors) == 1
    assert "step 2" in result.errors[0] and "Files in this Session" in result.errors[0]


@pytest.mark.parametrize(
    ("script", "reason"),
    [
        ("const track = content?.children?.[0];", "by position"),
        ("const row = card.children[1];", "by position"),
        ("async () => { await window.__GBDRAW_APP__.refreshLinearRecordSelectors(); }", "reads the records"),
        ("document.querySelector('.settings-pane').style.visibility = 'hidden';", "Settings pane"),
    ],
)
def test_capture_scripts_find_elements_by_identity_and_leave_the_record_read_to_the_tool(
    script: str, reason: str
) -> None:
    """A capture finds a control by its identity, not by its place among its siblings, which a new
    control shifts (OV-3366: the Layout crops showed Figure size instead of Track Layout). It does not
    read the records itself, nor click Inspect source records: the tool reads them after the replay
    for every capture (OV-3139)."""

    def step(action: dict[str, str]) -> dict[str, object]:
        return {"title": "Layout", "body": "Open Layout.", "operations": [{"capture": {"actions": [action]}}]}

    result = ValidationResult(errors=[], warnings=[])
    capture_tool.add_step_text_validation(result, _manual_context(2), step({"type": "evaluate", "script": script}))
    capture_tool.add_step_text_validation(
        result, _manual_context(3), step({"type": "click", "selector": "button:text-is('Inspect source records')"})
    )
    capture_tool.add_step_text_validation(
        result,
        _manual_context(4),
        step({"type": "evaluate", "script": "document.querySelector('[data-layout-options]').dataset.capture = 'x';"}),
    )

    assert [error.split(":")[0] for error in result.errors] == [_manual_context(2).label, _manual_context(3).label]
    assert reason in result.errors[0]


def test_capture_finds_its_card_by_section_identity_not_by_text() -> None:
    """A capture names its card by the section key, not by text: "Layout" matched the first card
    that says "layout" anywhere, so majanivirus step 4 showed Advanced comparison and layout in
    place of the Layout card its caption describes."""

    def step(capture: dict[str, object]) -> dict[str, object]:
        return {"title": "Layout", "body": "Open Layout.", "operations": [{"capture": capture}]}

    result = ValidationResult(errors=[], warnings=[])
    capture_tool.add_step_text_validation(result, _manual_context(2), step({"cardTitle": "Layout"}))
    capture_tool.add_step_text_validation(
        result, _manual_context(3), step({"selector": ".card[data-settings-section='layout']"})
    )
    capture_tool.add_step_text_validation(
        result, _manual_context(4), step({"openSelect": {"selectedLabel": "orchid", "cardTitle": "Colors", "selectIndex": 1}})
    )
    capture_tool.add_step_text_validation(
        result, _manual_context(5), step({"openSelect": [{"selectedLabel": "orchid", "selector": "select[aria-label='Palette']"}]})
    )

    assert [error.split(":")[0] for error in result.errors] == [_manual_context(2).label, _manual_context(4).label]
    assert "cardTitle" in result.errors[0] and "selectIndex" in result.errors[1]


def _text_box(x: float, y: float, width: float, height: float, text: str = "text") -> dict[str, object]:
    return {"x": x, "y": y, "width": width, "height": height, "text": text}


_EDGE_PAGE = {"width": 2000.0, "height": 2000.0}
_EDGE_PADDING = {"top": 20.0, "right": 20.0, "bottom": 20.0, "left": 20.0}
# The crop of a subject at x 520..780, y 520..680 with 20 px of padding, far from the page edge.
_EDGE_CLIP = {"x": 500.0, "y": 500.0, "width": 300.0, "height": 200.0}


def test_crop_edge_moves_to_the_nearest_position_that_cuts_no_text() -> None:
    """OV-3410: a crop edge that crosses a line of text or a control moves to the nearer position
    that crosses none: outward within CROP_EDGE_REACH, or inward within the padding, never into the
    subject; outward on a tie."""

    margin = capture_tool.CROP_EDGE_MARGIN
    reach = capture_tool.CROP_EDGE_REACH

    def cleared(*boxes: dict[str, object]) -> dict[str, float]:
        return capture_tool.clear_crop_edges(_EDGE_CLIP, list(boxes), _EDGE_PADDING, _EDGE_PAGE)

    # A heading mostly outside the top edge: the edge moves up over it.
    assert cleared(_text_box(510, 494, 120, 18, "DEFAULT COLORS")) == {
        "x": 500.0, "y": 494 - margin, "width": 300.0, "height": 206 + margin,
    }
    # A heading mostly inside the padding: the edge moves down under it, the subject stays whole.
    assert cleared(_text_box(510, 490, 120, 18, "DEFAULT COLORS"))["y"] == 508 + margin
    # Equally near both ways: outward.
    assert cleared(_text_box(510, 492, 120, 16, "Reset"))["y"] == 492 - margin
    # A glyph 1 px over the edge is a cut too (the descenders of a label above a field).
    assert cleared(_text_box(510, 483, 120, 18, "Max hits per protein"))["y"] == 501 + margin
    # A row that reaches into the subject: the edge cannot move up into it, so it moves down.
    assert cleared(_text_box(510, 676, 200, 30, "Ruler on Axis"))["height"] == 206 + margin
    # Text beside the crop, and text inside it, are not cut.
    assert cleared(_text_box(820, 490, 100, 18), _text_box(530, 530, 100, 18)) == _EDGE_CLIP
    # Text taller than the reach: the edge moves in, under it, into the padding only.
    assert cleared(_text_box(510, 508 - reach - 10, 120, reach + 10))["y"] == 508 + margin
    # Labels that touch one another: the edge moves past each, as far as the first gap.
    ticks = [_text_box(600, 676 + 12 * index, 20, 14, "3 kbp") for index in range(3)]
    assert cleared(*ticks)["height"] == 214 + margin
    # Moving the left edge out widens the top edge, which then meets text it did not cross before.
    label = _text_box(480, 600, 40, 16, "kbp")
    title = _text_box(460, 492, 30, 16, "Homo")
    widened = cleared(label, title)
    assert widened["x"] == 460 - margin and widened["y"] == 492 - margin
    assert capture_tool.crop_edge_cuts(widened, [label, title], _EDGE_PAGE) == []
    # The page edge cuts nothing: a crop at the page edge stays there.
    at_edge = {"x": 0.0, "y": 0.0, "width": 300.0, "height": 200.0}
    assert capture_tool.clear_crop_edges(at_edge, [_text_box(-5, -5, 50, 20)], _EDGE_PADDING, _EDGE_PAGE) == at_edge


def test_crop_guard_refuses_an_edge_that_still_crosses_text() -> None:
    """OV-3410: text the edge can neither take in (beyond the reach) nor leave out (it reaches into
    the subject) leaves the crop as it was, and the capture refuses it, naming the text."""

    reach = capture_tool.CROP_EDGE_REACH
    blocked = _text_box(510, 500 - reach - 10, 120, reach + 40, "Homo sapiens")
    cleared = capture_tool.clear_crop_edges(_EDGE_CLIP, [blocked], _EDGE_PADDING, _EDGE_PAGE)
    assert cleared == _EDGE_CLIP
    with pytest.raises(RuntimeError, match="top edge cuts .*Homo sapiens"):
        capture_tool.assert_crop_edges_cut_no_text(cleared, [blocked], _EDGE_PAGE)
    capture_tool.assert_crop_edges_cut_no_text(_EDGE_CLIP, [_text_box(530, 530, 100, 18)], _EDGE_PAGE)


def test_capture_parks_the_pointer_then_moves_crop_edges_off_text(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-3411: the pointer is parked after the capture's actions, so no crop shows a hover the last
    click left (BGC 7-2's Skip row); OV-3410: the screenshot takes the crop whose edges cut no text."""

    events: list[object] = []

    class Mouse:
        def move(self, x: float, y: float) -> None:
            events.append("park pointer")

    class FakePage:
        mouse = Mouse()

        def evaluate(self, script: str, *args: object) -> list[object]:
            return []

        def screenshot(self, **kwargs: object) -> bytes:
            events.append(("screenshot", kwargs["clip"]))
            buffer = io.BytesIO()
            Image.new("RGB", (2, 2)).save(buffer, "PNG")
            return buffer.getvalue()

    heading = _text_box(10, 95, 80, 10, "Palette")
    monkeypatch.setattr(capture_tool, "prepare_capture_page", lambda *args: events.append("prepare"))
    monkeypatch.setattr(capture_tool, "page_scroll_dimensions", lambda page: dict(_EDGE_PAGE))
    # raising=False: on a tool without the measurement the test fails on what it checks.
    monkeypatch.setattr(
        capture_tool,
        "visible_text_boxes",
        lambda page, region: events.append("measure text") or [heading],
        raising=False,
    )
    monkeypatch.setattr(capture_tool, "output_path_for_operation", lambda operation: tmp_path / "crop.webp")
    operation = {"capture": {"boundingBox": {"x": 0, "y": 100, "width": 200, "height": 100}}}

    capture_tool.capture_operation(FakePage(), {}, operation, "http://127.0.0.1", 90)

    assert events[:3] == ["prepare", "park pointer", "measure text"]
    assert events[-1] == (
        "screenshot",
        {"x": 0.0, "y": 95 - 4.0, "width": 200.0, "height": 105 + 4.0},
    )


def test_highlight_frame_moves_in_where_text_sits_closer_than_its_reach() -> None:
    """A highlight frame reaches HIGHLIGHT_FRAME_OFFSET + HIGHLIGHT_FRAME_WIDTH px outside its
    element, less where a box of text outside the element sits closer, and never covers it."""

    field = {"x": 100.0, "y": 100.0, "width": 120.0, "height": 24.0}
    room = capture_tool.highlight_frame_room
    style = capture_tool.highlight_frame_style
    width = capture_tool.HIGHLIGHT_FRAME_WIDTH
    clearance = capture_tool.HIGHLIGHT_FRAME_CLEARANCE

    assert room(field, []) == float("inf")
    assert style(float("inf")) == {
        "offset": capture_tool.HIGHLIGHT_FRAME_OFFSET,
        "glow": capture_tool.HIGHLIGHT_FRAME_GLOW,
    }
    # The field's own control and value overlap it: they do not limit the frame.
    own = [{"x": 100.0, "y": 100.0, "width": 120.0, "height": 24.0, "text": "input: near"},
           {"x": 108.0, "y": 104.0, "width": 30.0, "height": 16.0, "text": "0.4"}]
    assert room(field, own) == float("inf")
    # A label whose descenders end 3 px above the field: the frame stops 1 px short of them.
    label = {"x": 100.0, "y": 79.0, "width": 140.0, "height": 18.0, "text": "Head Length Ratio"}
    assert room(field, [label]) == 3.0
    near = style(3.0)
    assert near["offset"] + width == 3.0 - clearance and near["glow"] == 3.0 - clearance
    # Text touching the field: the frame lies on the field's own edge, with no glow.
    assert style(0.0) == {"offset": -width, "glow": 0.0}
    # A box off a corner limits the frame by the larger of its two gaps.
    corner = {"x": 230.0, "y": 60.0, "width": 40.0, "height": 20.0, "text": "x"}
    assert room(field, [corner]) == 20.0
