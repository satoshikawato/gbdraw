"""Shared visible-UI mechanics for five whole linear BGC records."""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import (
    assert_gui_bgc_collinear_svg,
    assert_gui_bgc_similarity_groups_svg,
    assert_plain_gui_bgc_svg,
    inspect_gui_bgc_losatp_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    GUI_BGC_FIXTURES,
)
from flows.web_capture import (
    assert_fixture_identity,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    open_linear_comparison_disclosure,
    select_linear_losat_mode,
    set_feature_search_visible,
    wait_for_worker,
)


BGC_DEFAULT_COLORS = GUI_BGC_FIXTURES[0][0].with_name(
    "BGC0000708-BGC0000713_default_colors.tsv"
)
BGC_SPECIFIC_COLORS = GUI_BGC_FIXTURES[0][0].with_name(
    "BGC0000708-BGC0000713_specific_colors.tsv"
)
BGC_GENE_PRIORITY = (
    GUI_BGC_FIXTURES[0][0].parents[1]
    / "shared"
    / "cds_gene_qualifier_priority.tsv"
)
BGC_RECORD_PRESENTATION = (
    (
        "<i>Streptomyces lividus</i> CBS 844.73",
        "Lividomycin biosynthetic gene cluster",
        False,
    ),
    (
        "<i>Streptomyces fradiae</i> ATCC 10745",
        "Neomycin biosynthetic gene cluster",
        False,
    ),
    (
        "<i>Streptomyces fradiae</i> MCIMB 8233",
        "Neomycin biosynthetic gene cluster",
        False,
    ),
    (
        "<i>Streptomyces rimosus</i> subsp. <i>paromomycinus</i> NRRL 2455",
        "Paromomycin biosynthetic gene cluster",
        False,
    ),
    (
        "<i>Streptomyces ribosidificus</i> ATCC 21294",
        "Ribostamycin biosynthetic gene",
        True,
    ),
)
EXPECTED_MEMBER_PROTEIN_IDS = (
    "CAG38695.1",
    "CAF33310.1",
    "CAH58688.1",
    "CAF32372.1",
    "CAG34720.1",
)


@dataclass(frozen=True)
class BgcLosatpResult:
    """Artifacts and semantic reports from one complete BGC journey."""

    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    tsv_download: dict[str, Any] | None
    fasta_download: dict[str, Any] | None
    popup: dict[str, Any]
    group_count: int


def assert_bgc_fixtures() -> None:
    """Require five unchanged, whole, naturally linear BGC records."""

    for path, size, digest, record_id, length, _organism in GUI_BGC_FIXTURES:
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        records = list(SeqIO.parse(path, "genbank"))
        if len(records) != 1:
            raise AssertionError(f"Expected one complete record in {path.name}")
        record = records[0]
        if record.id != record_id or len(record) != length:
            raise AssertionError(
                f"Unexpected BGC source identity: {record.id} ({len(record)} bp)"
            )
        if str(record.annotations.get("topology", "")).lower() != "linear":
            raise AssertionError(f"{record_id} is not annotated as a linear record")
        if len([feature for feature in record.features if feature.type == "CDS"]) == 0:
            raise AssertionError(f"{record_id} has no annotated CDS features")


def _set_bgc_inputs(page: Page) -> None:
    linear = page.get_by_role("button", name="Linear", exact=True)
    linear.click()
    expect(linear).to_have_attribute("aria-pressed", "true")
    page.get_by_role("radio", name="GenBank", exact=True).check()
    expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
        "Current: No comparison"
    )

    add_sequence = page.get_by_role("button", name="Add sequence", exact=True)
    expect(add_sequence).to_have_count(2)
    for _ in range(4):
        add_sequence.first.click()
    for index, fixture in enumerate(GUI_BGC_FIXTURES, start=1):
        page.get_by_test_id(f"linear-genbank-{index}").set_input_files(fixture[0])

    selected_files = page.get_by_role(
        "group", name="GenBank File selection", exact=True
    )
    expect(selected_files).to_have_count(5)
    for index, fixture in enumerate(GUI_BGC_FIXTURES):
        expect(selected_files.nth(index)).to_contain_text(fixture[0].name)

    for index, (definition, subtitle, reverse) in enumerate(
        BGC_RECORD_PRESENTATION,
        start=1,
    ):
        record_options = page.get_by_role(
            "button",
            name=f"Record options for sequence {index}",
            exact=True,
        )
        record_options.click()
        definition_input = page.get_by_label(
            f"Definition for sequence {index}",
            exact=True,
        )
        definition_input.fill(definition)
        expect(definition_input).to_have_value(definition)
        subtitle_input = page.get_by_label(
            f"Subtitle / title for sequence {index}",
            exact=True,
        )
        subtitle_input.fill(subtitle)
        expect(subtitle_input).to_have_value(subtitle)
        reverse_control = page.get_by_label(
            f"Reverse complement for sequence {index}",
            exact=True,
        )
        if reverse:
            reverse_control.check()
            expect(reverse_control).to_be_checked()
        else:
            expect(reverse_control).not_to_be_checked()
        record_options.click()


def _set_gallery_quality_presentation(
    page: Page,
    *,
    title: str,
    match_gallery_definitions: bool,
) -> None:
    colors = page.get_by_label("Colors", exact=True)
    colors.click()
    palette = page.get_by_label("Palette", exact=True)
    palette.select_option("orange")
    expect(palette).to_have_value("orange")
    page.get_by_label("Override File (-d)", exact=True).set_input_files(
        BGC_DEFAULT_COLORS
    )
    page.get_by_label("Specific Table (-t)", exact=True).set_input_files(
        BGC_SPECIFIC_COLORS
    )
    colors.click()

    labels = page.get_by_label("Labels", exact=True)
    labels.click()
    show_labels = page.get_by_label("Show Labels", exact=True)
    show_labels.select_option("first")
    expect(show_labels).to_have_value("first")
    page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
        BGC_GENE_PRIORITY
    )
    label_font_size = page.get_by_label("Label Font Size", exact=True)
    label_font_size.fill("18")
    expect(label_font_size).to_have_value("18")
    label_placement = page.get_by_label("Label Placement", exact=True)
    label_placement.select_option("above_feature")
    expect(label_placement).to_have_value("above_feature")
    label_rotation = page.get_by_label("Label Rotation", exact=True)
    label_rotation.fill("45")
    expect(label_rotation).to_have_value("45")
    labels.click()

    features = page.get_by_label("Features", exact=True)
    features.click()
    feature_height = page.get_by_label("Feature Height", exact=True)
    feature_height.fill("75")
    expect(feature_height).to_have_value("75")
    block_color_mode = page.get_by_label("Block stroke color mode", exact=True)
    block_color_mode.select_option("color")
    block_stroke_color = page.get_by_label("Block stroke color", exact=True)
    block_stroke_color.fill("#262626")
    expect(block_stroke_color).to_have_value("#262626")
    block_stroke_width = page.get_by_label("Block Stroke Width", exact=True)
    block_stroke_width.fill("2")
    expect(block_stroke_width).to_have_value("2")
    line_stroke_width = page.get_by_label("Line Stroke Width", exact=True)
    line_stroke_width.fill("2")
    expect(line_stroke_width).to_have_value("2")
    features.click()

    axis_and_scale = page.get_by_label("Axis & Scale", exact=True)
    axis_and_scale.click()
    show_scale = page.get_by_label("Show Coordinate Scale (Linear)", exact=True)
    show_scale.check()
    expect(show_scale).to_be_checked()
    scale_style = page.get_by_label("Linear scale style", exact=True)
    scale_style.select_option("ruler")
    expect(scale_style).to_have_value("ruler")
    axis_stroke_width = page.get_by_label("Axis Stroke Width", exact=True)
    axis_stroke_width.fill("5")
    expect(axis_stroke_width).to_have_value("5")
    axis_and_scale.click()

    title_and_legend = page.get_by_label("Title & Legend", exact=True)
    title_and_legend.click()
    page.get_by_label("Plot Title", exact=True).fill(title)
    page.get_by_label("Plot Title Position", exact=True).select_option("bottom")
    page.get_by_label("Legend Position", exact=True).select_option("bottom")

    if match_gallery_definitions:
        definition_styles = page.get_by_text(
            "Definition Line Styles", exact=True
        ).locator("xpath=../following-sibling::div[1]")

        def set_definition_line(
            label: str,
            *,
            size: str,
            bold: bool = False,
            color: str | None = None,
        ) -> None:
            line_label = definition_styles.get_by_text(label, exact=True)
            size_input = line_label.locator("xpath=following-sibling::input[1]")
            size_input.fill(size)
            expect(size_input).to_have_value(size)
            weight_controls = size_input.locator("xpath=following-sibling::div[1]")
            weight_controls.get_by_role(
                "button", name="Bold" if bold else "Normal", exact=True
            ).click()
            if color is not None:
                fill_input = page.get_by_label(
                    f"{label} definition line fill value", exact=True
                )
                fill_input.fill(color)
                expect(fill_input).to_have_value(color)

        set_definition_line("Name / Species", size="20", bold=True)
        set_definition_line("Subtitle", size="20")
        set_definition_line("Accession", size="20", color="#7b7c7d")
        set_definition_line("Length / Coord.", size="20", color="#7b7c7d")
    title_and_legend.click()

    track_layout = page.get_by_label("Track Layout", exact=True)
    track_layout.select_option("middle")
    if match_gallery_definitions:
        lock_definition = page.get_by_label("Lock Definition Column", exact=True)
        lock_definition.check()
        expect(lock_definition).to_be_checked()


def _configure_losatp(
    page: Page, *, mode: str, output_prefix: str
) -> tuple[Any, Any]:
    commands = page.get_by_role(
        "group", name="Set all adjacent comparisons", exact=True
    )
    run_losat = commands.get_by_role(
        "button", name="Run LOSAT for all adjacent pairs", exact=True
    )
    run_losat.click()
    expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
        "Current: Run LOSAT for all adjacent pairs"
    )

    settings = open_linear_comparison_disclosure(
        page,
        "settings",
        "Comparison Settings",
    )
    select_linear_losat_mode(
        settings,
        label="LOSATP",
        mode_key="blastp",
    )
    losatp_mode = settings.get_by_role(
        "combobox", name="LOSATP mode", exact=True
    )
    losatp_mode.select_option("pairwise")
    expect(losatp_mode).to_have_value("pairwise")
    match_style = settings.get_by_label("Pairwise Match Style", exact=True)
    match_style.select_option("curve")
    expect(match_style).to_have_value("curve")
    losatp_mode.select_option(mode)
    expect(losatp_mode).to_have_value(mode)

    settings.get_by_label(
        "Linear comparison minimum bitscore", exact=True
    ).fill("50")
    settings.get_by_label(
        "Linear comparison maximum e-value", exact=True
    ).fill("0.01")
    settings.get_by_label(
        "Linear comparison minimum identity", exact=True
    ).fill("30")
    settings.get_by_label(
        "Linear comparison minimum alignment length", exact=True
    ).fill("0")

    if mode == "collinear":
        settings.get_by_label("Collinear max unit gap", exact=True).fill("1")
        settings.get_by_label(
            "Collinear minimum block genes", exact=True
        ).fill("2")
        settings.get_by_label("Collinear color mode", exact=True).select_option(
            "orientation_identity"
        )
        settings.get_by_label(
            "Collinear evidence scope", exact=True
        ).select_option("adjacent")

    advanced = open_linear_comparison_disclosure(
        page,
        "advanced",
        "Advanced comparison and layout",
    )
    execution = advanced.get_by_role(
        "combobox", name="LOSAT execution", exact=True
    )
    execution.select_option("serial")
    total_threads = advanced.get_by_role(
        "combobox", name="LOSAT total threads", exact=True
    )
    total_threads.select_option("1")
    parallel_runs = advanced.get_by_role(
        "combobox", name="LOSAT parallel runs", exact=True
    )
    parallel_runs.select_option("1")
    threads_per_run = advanced.get_by_role(
        "combobox", name="LOSAT threads per run", exact=True
    )
    if threads_per_run.is_enabled():
        threads_per_run.select_option("1")

    if mode == "collinear":
        advanced.get_by_label("Collinear diagonal drift", exact=True).fill("1")
        advanced.get_by_label("Collinear merge conflicts", exact=True).fill("0")

    page.get_by_label("Output Prefix", exact=True).fill(output_prefix)
    return settings, advanced


def _inspect_popup(page: Page, *, mode: str) -> dict[str, Any]:
    match_kind = "orthogroup" if mode == "orthogroup" else "collinear"
    first_match = page.get_by_role(
        "region", name="Result Preview", exact=True
    ).locator(f'[data-match-kind="{match_kind}"]').first
    expect(first_match).to_be_visible()
    expect(first_match).to_have_attribute("role", "button")
    first_match.click()
    popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
    expect(popup).to_be_visible()
    text = popup.inner_text()
    normalized_text = text.casefold()
    required = (
        ("Similarity group ID", "Members")
        if mode == "orthogroup"
        else ("Block ID", "Anchors", "Orientation")
    )
    for token in required:
        if token.casefold() not in normalized_text:
            raise AssertionError(
                f"LOSATP popup is missing {token!r}: {text!r}"
            )
    return {"text": text, "mode": mode}


def _open_orthogroup_alignment_target(page: Page, orthogroup_id: str) -> Any:
    target = page.evaluate(
        """
        (orthogroupId) => {
          const app = window.__GBDRAW_APP__;
          const svg = document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement
            || document.querySelector('svg');
          const features = Array.isArray(app?.extractedFeatures)
            ? app.extractedFeatures
            : [];
          const feature = features.find(
            (item) => String(item?.orthogroupId || '').trim() === orthogroupId
          );
          const featureId = String(feature?.svg_id || '').trim();
          const element = Array.from(
            svg?.querySelectorAll('[data-gbdraw-feature-id]') || []
          ).find(
            (candidate) => String(
              candidate.getAttribute('data-gbdraw-feature-id') || ''
            ).trim() === featureId
          );
          if (!element) throw new Error(`Rendered ${orthogroupId} feature was not found`);
          element.scrollIntoView({ block: 'center', inline: 'center' });
          const rect = element.getBoundingClientRect();
          return {
            featureId,
            x: rect.left + rect.width / 2,
            y: rect.top + rect.height / 2
          };
        }
        """,
        orthogroup_id,
    )
    page.evaluate(
        """
        ({ featureId, x, y }) => {
          const svg = document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement
            || document.querySelector('svg');
          const element = Array.from(
            svg?.querySelectorAll('path[data-gbdraw-feature-id], polygon[data-gbdraw-feature-id], rect[data-gbdraw-feature-id]') || []
          ).find(
            (candidate) => String(
              candidate.getAttribute('data-gbdraw-feature-id') || ''
            ).trim() === featureId
          );
          if (!element) throw new Error(`Rendered feature ${featureId} was not found`);
          element.dispatchEvent(new MouseEvent('click', {
            bubbles: true,
            cancelable: true,
            clientX: x,
            clientY: y,
            view: window
          }));
        }
        """,
        target,
    )
    page.wait_for_function(
        """
        ({ featureId, orthogroupId }) => {
          const feature = window.__GBDRAW_APP__?.clickedFeature;
          return feature?.svg_id === featureId
            && feature?.orthogroupId === orthogroupId;
        }
        """,
        arg={"featureId": target["featureId"], "orthogroupId": orthogroup_id},
    )
    popup = page.locator(".feature-popup")
    expect(popup).to_be_visible()
    expect(popup).to_contain_text(orthogroup_id)
    expect(popup.get_by_role("button", name=re.compile(r"Align"))).to_be_visible()
    return popup


def _align_to_orthogroup(page: Page, popup: Any, orthogroup_id: str) -> None:
    previous_run = page.evaluate(
        "() => String(window.__GBDRAW_APP__?.lastRunInfo?.startedAtIso || '')"
    )
    popup.get_by_role("button", name=re.compile(r"Align")).click()
    page.wait_for_function(
        """
        ({ orthogroupId, previousRun }) => {
          const app = window.__GBDRAW_APP__;
          return !app?.processing
            && app?.selectedOrthogroupAlignmentFeature === orthogroupId
            && String(app?.lastRunInfo?.startedAtIso || '') !== previousRun;
        }
        """,
        arg={"orthogroupId": orthogroup_id, "previousRun": previous_run},
        timeout=ACTION_TIMEOUT_MS,
    )


def _download_text_button(
    page: Page,
    *,
    button: Any,
    download_dir: Path,
    expected_name: str | None = None,
) -> tuple[Path, str]:
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        button.click()
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"Browser download failed: {download.failure()}")
    suggested = download.suggested_filename
    if expected_name is not None and suggested != expected_name:
        raise AssertionError(
            f"Expected download {expected_name}, found {suggested}"
        )
    target = download_dir / (expected_name or suggested)
    download.save_as(target)
    return target, suggested


def _validate_raw_losatp(path: Path, *, expected_name: str) -> dict[str, Any]:
    if path.name != expected_name:
        raise AssertionError(f"Unexpected LOSATP TSV name: {path.name}")
    rows = [
        line.split("\t")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line and not line.startswith("#")
    ]
    if not rows or any(len(row) != 12 for row in rows):
        raise AssertionError("LOSATP output must contain non-empty 12-column rows")
    return {"filename": path.name, "bytes": path.stat().st_size, "rows": len(rows)}


def _expected_member_translations() -> dict[str, str]:
    expected_ids = set(EXPECTED_MEMBER_PROTEIN_IDS)
    translations: dict[str, str] = {}
    for fixture, *_ in GUI_BGC_FIXTURES:
        record = next(SeqIO.parse(fixture, "genbank"))
        for feature in record.features:
            if feature.type != "CDS":
                continue
            protein_ids = feature.qualifiers.get("protein_id", [])
            if len(protein_ids) != 1 or protein_ids[0] not in expected_ids:
                continue
            protein_id = protein_ids[0]
            translation = "".join(feature.qualifiers.get("translation", []))
            if not translation or protein_id in translations:
                raise AssertionError(
                    f"Ambiguous source translation for {protein_id} in {fixture.name}"
                )
            translations[protein_id] = translation
    if set(translations) != expected_ids:
        raise AssertionError(
            "Raw BGC records do not contain the expected member translations: "
            f"{sorted(translations)!r}"
        )
    return translations


def _validate_fasta(path: Path, *, suggested_name: str) -> dict[str, Any]:
    records = list(SeqIO.parse(path, "fasta"))
    expected = _expected_member_translations()
    downloaded = {record.id: str(record.seq) for record in records}
    if len(downloaded) != len(records) or set(downloaded) != set(expected):
        raise AssertionError(
            "LOSATP member FASTA changed: "
            f"{sorted(downloaded)!r} ({suggested_name})"
        )
    mismatches = [
        protein_id
        for protein_id, translation in expected.items()
        if downloaded[protein_id] != translation
    ]
    if mismatches:
        raise AssertionError(
            "LOSATP member FASTA sequences differ from the raw GenBank CDS "
            f"translations: {mismatches!r} ({suggested_name})"
        )
    verified_members = [
        {
            "proteinId": protein_id,
            "length": len(expected[protein_id]),
            "sha256": hashlib.sha256(
                expected[protein_id].encode("ascii")
            ).hexdigest(),
        }
        for protein_id in EXPECTED_MEMBER_PROTEIN_IDS
    ]
    return {
        "filename": path.name,
        "suggestedFilename": suggested_name,
        "bytes": path.stat().st_size,
        "records": len(records),
        "verifiedMembers": verified_members,
    }


def capture_bgc_losatp(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
    *,
    scenario_id: str,
    mode: str,
    output_prefix: str,
    raw_tsv_name: str | None,
    screenshot_names: Mapping[str, str],
    include_plain_result: bool,
    download_member_fasta: bool,
    separate_strands: bool | None = None,
    assert_final: Callable[[dict[str, Any]], None] | None = None,
    align_orthogroup_id: str | None = None,
    match_gallery_definitions: bool = False,
) -> BgcLosatpResult:
    """Run one real serial LOSATP journey from five complete source files."""

    assert_bgc_fixtures()
    if mode not in {"orthogroup", "collinear"}:
        raise ValueError(f"Unsupported LOSATP documentation mode: {mode}")
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    tsv_report: dict[str, Any] | None = None
    fasta_report: dict[str, Any] | None = None

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _set_bgc_inputs(page)
        if separate_strands is not None:
            separate_strands_control = page.get_by_label(
                "Separate Strands", exact=True
            )
            if separate_strands:
                separate_strands_control.check()
                expect(separate_strands_control).to_be_checked()
            else:
                separate_strands_control.uncheck()
                expect(separate_strands_control).not_to_be_checked()
        if "input" in screenshot_names:
            name = screenshot_names["input"]
            fifth_record = page.get_by_role(
                "group", name="Linear sequence 5", exact=True
            )
            fifth_file = fifth_record.get_by_role(
                "group", name="GenBank File selection", exact=True
            )
            reverse_complement = page.get_by_label(
                "Reverse complement for sequence 5", exact=True
            )
            fifth_options = page.get_by_role(
                "button",
                name="Record options for sequence 5",
                exact=True,
            )
            fifth_options.click()
            reverse_complement.evaluate(
                "(element) => element.scrollIntoView({ block: 'center' })"
            )
            expect(fifth_file).to_contain_text("BGC0000713.gbk")
            expect(fifth_file).to_be_in_viewport()
            expect(reverse_complement).to_be_checked()
            expect(reverse_complement).to_be_in_viewport()
            screenshot_bytes[name] = capture_screenshot(
                page, output_paths[name], "Linear"
            )
            fifth_options.click()

        if include_plain_result:
            generate_and_inspect(
                page, inspect_gui_bgc_losatp_svg, assert_plain_gui_bgc_svg
            )
            fit_complete_linear_preview(page, target_zoom="30%")
            name = screenshot_names["plain"]
            screenshot_bytes[name] = capture_screenshot(
                page, output_paths[name], "Linear"
            )

        title = (
            "LOSATP Similarity groups across five whole BGC records"
            if mode == "orthogroup"
            else "LOSATP Collinear blocks across five whole BGC records"
        )
        _set_gallery_quality_presentation(
            page,
            title=title,
            match_gallery_definitions=match_gallery_definitions,
        )
        settings, advanced = _configure_losatp(
            page,
            mode=mode,
            output_prefix=output_prefix,
        )

        validator = assert_final or (
            assert_gui_bgc_similarity_groups_svg
            if mode == "orthogroup"
            else assert_gui_bgc_collinear_svg
        )
        final_report = generate_and_inspect(
            page, inspect_gui_bgc_losatp_svg, validator
        )
        if mode == "orthogroup":
            group_count = int(
                page.evaluate(
                    """
                    () => {
                      const app = window.__GBDRAW_APP__;
                      if (Array.isArray(app?.orthogroups)) return app.orthogroups.length;
                      return Number(app?.orthogroupCount?.value ?? app?.orthogroupCount ?? 0);
                    }
                    """
                )
            )
            if group_count != 23:
                group_debug = page.evaluate(
                    """
                    () => {
                      const app = window.__GBDRAW_APP__;
                      return {
                        mode: app?.losat?.blastp?.mode,
                        program: app?.losatProgram,
                        orthogroupsIsArray: Array.isArray(app?.orthogroups),
                        orthogroupsLength: Array.isArray(app?.orthogroups)
                          ? app.orthogroups.length
                          : null,
                        orthogroupCount: app?.orthogroupCount,
                        comparisonMode: app?.linearComparisonMode,
                        lastRunInfo: app?.lastRunInfo
                      };
                    }
                    """
                )
                raise AssertionError(
                    "Expected 23 stable LOSATP groups, "
                    f"found {group_count}: {group_debug!r}"
                )
            telemetry = page.evaluate(
                "() => window.__GBDRAW_APP__?.lastRunInfo?.losatTelemetry || {}"
            )
            expected_jobs = len(GUI_BGC_FIXTURES) ** 2
            if (
                telemetry.get("totalPairs") != expected_jobs
                or telemetry.get("uniqueJobs") != expected_jobs
            ):
                raise AssertionError(
                    "Similarity groups did not use all-vs-all LOSATP evidence: "
                    f"{telemetry!r}"
                )
        else:
            group_count = len(
                {
                    str(match.get("blockId") or "")
                    for match in final_report["comparisonMatches"]
                    if match.get("blockId")
                }
            )

        set_feature_search_visible(page, visible=False)
        fit_complete_linear_preview(page, target_zoom="40%")
        name = screenshot_names["settings"]
        losat_mode_group = settings.get_by_role(
            "group", name="LOSAT Mode", exact=True
        )
        losat_mode_group.evaluate(
            """
            (element) => {
              element.scrollIntoView({ block: 'center' });
              let owner = element.parentElement;
              while (owner) {
                owner.scrollLeft = 0;
                owner = owner.parentElement;
              }
              window.scrollTo(0, window.scrollY);
            }
            """
        )
        expect(losat_mode_group).to_be_in_viewport()
        screenshot_bytes[name] = capture_screenshot(
            page, output_paths[name], "Linear"
        )
        if align_orthogroup_id is not None:
            alignment_popup = _open_orthogroup_alignment_target(
                page, align_orthogroup_id
            )
            alignment_name = screenshot_names["align"]
            screenshot_bytes[alignment_name] = capture_screenshot(
                page, output_paths[alignment_name], "Linear"
            )
            _align_to_orthogroup(page, alignment_popup, align_orthogroup_id)
            result_region = page.get_by_role(
                "region", name="Result Preview", exact=True
            )
            final_report = inspect_gui_bgc_losatp_svg(result_region)
            validator(final_report)
            aligned_state = page.evaluate(
                "() => window.__GBDRAW_APP__?.selectedOrthogroupAlignmentFeature"
            )
            if aligned_state != align_orthogroup_id:
                raise AssertionError(
                    f"Expected alignment to {align_orthogroup_id}, found {aligned_state!r}"
                )
            fit_complete_linear_preview(page, target_zoom="40%")
        else:
            fit_complete_linear_preview(page, target_zoom="40%")
        if scenario_id == "H-GUI-07":
            raw_result = advanced.locator('[data-linear-raw-result]').first
            raw_result.evaluate(
                "(element) => element.scrollIntoView({ block: 'center' })"
            )
            expect(raw_result).to_contain_text("Raw result ready")
        name = screenshot_names["result"]
        screenshot_bytes[name] = capture_screenshot(
            page, output_paths[name], "Linear"
        )

        if raw_tsv_name is not None:
            raw_name = advanced.get_by_role(
                "textbox", name="Raw LOSAT filename for #1 to #2", exact=True
            )
            raw_button = advanced.get_by_role(
                "button", name="Save Raw LOSAT TSV for #1 to #2", exact=True
            )
            expect(raw_button).to_be_enabled()
            raw_name.fill(raw_tsv_name)
            raw_name.press("Tab")
            raw_path, _ = _download_text_button(
                page,
                button=raw_button,
                download_dir=download_dir,
                expected_name=raw_tsv_name,
            )
            tsv_report = _validate_raw_losatp(raw_path, expected_name=raw_tsv_name)

        popup_report = _inspect_popup(page, mode=mode)
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        popup.hover()
        page.mouse.wheel(0, 420)
        name = screenshot_names["popup"]
        screenshot_bytes[name] = capture_screenshot(
            page, output_paths[name], "Linear"
        )

        if download_member_fasta:
            group_row = popup.get_by_role("row").filter(has_text="og_1").first
            expect(group_row).to_be_visible()
            group_row.click()
            fasta_button = page.get_by_title(
                "Download all member amino-acid FASTA", exact=True
            ).first
            expect(fasta_button).to_be_enabled()
            fasta_path, suggested = _download_text_button(
                page,
                button=fasta_button,
                download_dir=download_dir,
            )
            normalized_fasta = download_dir / "collinear_members.fasta"
            if fasta_path != normalized_fasta:
                fasta_path.replace(normalized_fasta)
            fasta_report = _validate_fasta(
                normalized_fasta,
                suggested_name=suggested,
            )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        expected_svg_name = f"{output_prefix}.svg"
        svg_path, _ = _download_text_button(
            page,
            button=svg_button,
            download_dir=download_dir,
            expected_name=expected_svg_name,
        )
        svg_bytes = svg_path.stat().st_size
        if svg_bytes < 10_000:
            raise AssertionError(f"Downloaded BGC SVG is too small: {svg_bytes}")
        download_report = {
            "filename": svg_path.name,
            "bytes": svg_bytes,
        }

        capture.assert_clean()
    finally:
        capture.close()

    return BgcLosatpResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        tsv_download=tsv_report,
        fasta_download=fasta_report,
        popup=popup_report,
        group_count=group_count,
    )
