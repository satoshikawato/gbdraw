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
    assert (
        ':aria-label="$attrs[\'data-input-aria-label\'] || label" '
        '@change="handleFile"' in index
    )
    assert (
        'role="group" :aria-label="($attrs[\'data-input-aria-label\'] || label) '
        '+ \' selection\'"' in index
    )

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
    assert 'role="group" aria-label="Set all adjacent comparisons"' in index
    for name in (
        "Set no comparison",
        "Run LOSAT for all adjacent pairs",
        "Use uploaded BLAST TSV for all adjacent pairs",
    ):
        assert f'aria-label="{name}"' in index
    command_group = index[
        index.index('aria-label="Set all adjacent comparisons"'):
        index.index(
            'role="status"',
            index.index('aria-label="Set all adjacent comparisons"'),
        )
    ]
    assert "aria-pressed" not in command_group
    assert 'role="status" aria-live="polite"' in index
    assert 'data-linear-comparison-disclosure="settings"' in index
    assert 'aria-label="Comparison Settings"' in index
    assert 'data-linear-comparison-losat-mode' in index
    assert 'aria-label="LOSAT Mode"' in index
    assert 'data-linear-comparison-losat-mode-option="modeOption.key"' in index
    assert ':aria-pressed="modeOption.active"' in index
    assert '@click="setLinearComparisonLosatMode(modeOption.key)"' in index
    assert 'id="linear-comparison-losat-mode"' not in index
    assert 'data-linear-comparison-losatp-mode' in index
    assert 'aria-label="LOSATP mode"' in index
    assert '@change="setLinearComparisonLosatpMode($event.target.value)"' in index
    assert "sectionKeys.settings.includes('losatp-mode')" in index
    assert 'Comparison search method' not in index
    assert 'role="alert" aria-label="Generation Error"' in index
    assert 'aria-label="Add sequence" @click="addLinearSeq()"' in index


def test_linear_input_comparison_and_generate_controls_follow_semantic_dom_order() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")

    assert 'data-linear-record-list' in index
    assert (
        'v-for="displayRow in linearComparisonTimeline.rows"' in index
    )
    assert ':data-linear-display-row="displayRow.row"' in index
    assert (
        'v-for="({ sequence: seq, index: idx }) in displayRow.records"'
        in index
    )
    assert ':key="seq.uid"' in index
    assert ':data-linear-record-card="seq.uid"' in index
    assert ':data-linear-record-uid="seq.uid"' in index
    assert ':data-linear-record-options="seq.uid"' in index
    assert 'data-linear-comparison-card' in index
    assert 'data-linear-comparison-timeline' in index
    assert 'data-linear-comparison-disclosure="selected-pairs"' in index
    assert (
        ':data-linear-comparison-boundary="`${displayRow.boundaryAfter.'
        'upperRow}->${displayRow.boundaryAfter.lowerRow}`"' in index
    )
    assert (
        '<fieldset\n'
        '                                        v-for="pair in '
        'displayRow.boundaryAfter.pairs"' in index
    )
    assert ':data-edge-key="pair.edgeKey"' in index
    assert ':data-linear-comparison-pair-row="pair.edgeKey"' in index
    assert ':data-linear-comparison-pair-upload="pair.edgeKey"' in index
    assert '@change="setLinearComparisonGapAction(pair.edgeKey, ' in index
    assert (
        '@update:model-value="setLinearComparisonCardFile(pair.edgeKey, '
        '$event)"' in index
    )
    assert (
        '@change="setResolvedLinearComparisonLosatFilename('
        'pair.edgeKey, $event.target.value)"' in index
    )
    assert (
        '@click="deactivateResolvedLinearComparisonLosatFilename('
        'pair.edgeKey)"' in index
    )

    header_add = index.index('aria-label="Add sequence"')
    record_list = index.index("data-linear-record-list")
    comparison = index.index("data-linear-comparison-card")
    timeline = index.index("data-linear-comparison-timeline")
    basic = index.index('<div class="card basic-settings">')
    generate = index.index('<div class="generate-bar')
    advanced = index.index("data-linear-advanced-comparison")
    assert header_add < record_list < comparison < timeline < basic < generate < advanced
    assert index.count('aria-label="Generate Diagram" @click="runAnalysis"') == 1

    first_record_uploader = index.index('label="GenBank File"', record_list)
    first_record_options = index.index(':data-linear-record-options="seq.uid"')
    assert first_record_uploader < first_record_options < comparison

    assert 'aria-label="Selected pairs"' in index
    assert 'aria-label="Advanced comparison and layout"' in index
    assert (
        ':data-linear-unplaced-draft="unplaced.draft.id"' in index
    )
    assert "linearComparisonTimeline.unplacedDrafts" in index
    assert "Zipped adjacent pairs" in index
    assert "All adjacent-row pairs (cross-product)" in index
    assert ">Omit all</button>" in index

    assert "linearAdjacentComparisonRows" not in index
    assert "linearResolvedLosatEdges" not in index
    assert ">Adjacent gaps</" not in index
    assert ">Selected pairs and retained drafts</" not in index
    raw_results = index.index(">Raw LOSAT results</div>")
    assert advanced < raw_results


def test_file_uploader_exposes_a_native_keyboard_trigger() -> None:
    index = WEB_INDEX.read_text(encoding="utf-8")

    template_start = index.index(
        '<script type="text/x-template" id="file-uploader-template">'
    )
    template = index[template_start : index.index("</script>", template_start)]
    assert 'role="button"' in template
    assert 'tabindex="0"' in template
    assert ':aria-label="`Choose ${$attrs[\'data-input-aria-label\'] || label}`"' in template
    assert '@keydown.enter.prevent="$refs.input.click()"' in template
    assert '@keydown.space.prevent="$refs.input.click()"' in template
    assert '.upload-zone:focus-visible' in index


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
    ):
        assert f'aria-label="{label}"' in index

    assert (
        ':aria-label="`Raw LOSAT filename for #${pair.queryIndex + 1} '
        'to #${pair.subjectIndex + 1}`"' in index
    )
    assert (
        ':aria-label="`Save Raw LOSAT TSV for #${pair.queryIndex + 1} '
        'to #${pair.subjectIndex + 1}`"' in index
    )
    assert (
        ':data-input-aria-label="`BLAST TSV for #${pair.queryIndex + 1} to '
        '#${pair.subjectIndex + 1}`"' in index
    )
    assert 'label="BLAST TSV"' in index
    assert ':data-linear-comparison-pair-upload="pair.edgeKey"' in index
    assert ':data-input-aria-label="`BLAST TSV for #' in index
    assert ':test-id="`linear-genbank-${idx + 1}`"' in index
    assert ':data-testid="testId || undefined"' in index
    assert "'testId'" in components
    assert 'role="dialog"' in index
    assert 'aria-label="Pairwise match details"' in index
    assert "element.setAttribute('role', 'button')" in svg_actions
    assert "element.setAttribute('tabindex', '0')" in svg_actions
    assert "`Pairwise match ${index + 1}`" in svg_actions
    assert "e.key !== 'Enter' && e.key !== ' '" in svg_actions
