from pathlib import Path


WEB_INDEX = (
    Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "index.html"
)


def test_circular_legend_selector_exposes_every_supported_corner() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")

    for position in ("upper_left", "upper_right", "lower_left", "lower_right"):
        assert (
            f'<option value="{position}" v-if="mode === \'circular\'">'
            in index
        )


def test_first_circular_tutorial_controls_have_stable_accessible_selectors() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")

    assert 'aria-label="Circular" :aria-pressed="mode === \'circular\'"' in index
    assert 'aria-label="Linear" :aria-pressed="mode === \'linear\'"' in index
    assert 'ref="input" :multiple="multiple"' in index
    assert ':aria-label="label" @change="handleFile"' in index
    assert 'role="group" :aria-label="label + \' selection\'"' in index

    for label, control_id in (
        ("Output Prefix", "output-prefix"),
        ("Species", "circular-species"),
        ("Track Preset", "circular-track-preset"),
        ("Legend Position", "legend-position"),
        ("Label Mode", "circular-label-mode"),
    ):
        assert f'<label for="{control_id}"' in index, label
        assert f'id="{control_id}"' in index, label

    for name in ("Separate Strands", "Hide GC Content", "Hide GC Skew"):
        assert f'aria-label="{name}"' in index

    assert '<summary aria-label="Title &amp; Legend">' in index
    assert '<summary aria-label="Labels">' in index
    assert 'aria-label="Generate Diagram" @click="runAnalysis"' in index
    assert 'role="region" aria-label="Result Preview"' in index
    assert 'aria-label="SVG" @click="downloadSVG"' in index
    assert 'cursor-help" aria-hidden="true"' in index


def test_first_linear_tutorial_controls_have_stable_accessible_selectors() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")

    for label, control_id in (
        ("Track Layout", "linear-track-layout"),
        ("Show Labels", "linear-show-labels"),
    ):
        assert f'<label for="{control_id}"' in index, label
        assert f'id="{control_id}"' in index, label

    assert '<summary aria-label="Layout">' in index
    assert '<summary aria-label="Axis &amp; Scale">' in index
    assert 'aria-label="Show Coordinate Scale (Linear)"' in index
    assert 'aria-label="Linear scale style"' in index
    assert 'name="linear-comparison-global"' in index
    assert "> No comparison</label>" in index
    assert 'role="alert" aria-label="Generation Error"' in index
    assert 'aria-label="Add sequence" @click="addLinearSeq()"' in index


def test_losatn_tutorial_controls_and_match_popup_are_accessible() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")
    svg_actions = (
        WEB_INDEX.parent
        / "js"
        / "app"
        / "feature-editor"
        / "svg-actions.js"
    ).read_text(encoding="utf-8")
    components = (WEB_INDEX.parent / "js" / "components.js").read_text(
        encoding="utf-8"
    )

    for label in (
        "LOSAT execution",
        "LOSAT total threads",
        "LOSAT parallel runs",
        "LOSAT threads per run",
        "LOSATN task",
        "Raw LOSAT filename",
        "Save Raw LOSAT TSV",
    ):
        assert f'aria-label="{label}"' in index

    assert ':test-id="`linear-genbank-${idx + 1}`"' in index
    assert ':data-testid="testId || undefined"' in index
    assert "'testId'" in components
    assert 'role="dialog"' in index
    assert 'aria-label="Pairwise match details"' in index
    assert "element.setAttribute('role', 'button')" in svg_actions
    assert "element.setAttribute('tabindex', '0')" in svg_actions
    assert "`Pairwise match ${index + 1}`" in svg_actions
    assert "e.key !== 'Enter' && e.key !== ' '" in svg_actions
