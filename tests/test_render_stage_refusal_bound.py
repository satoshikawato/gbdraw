"""The warm-equals-cold harness renders a label-fit refusal the same on any host (OV-3483).

A Circular refusal computes its rows until a wall-clock bound
(``SWITCH_SEARCH_SECONDS``); a loaded host reaches it sooner, so a warm and a
cold render of the same request could list different rows. The harness renders
without that bound, so the message it compares does not depend on the load.
"""

from __future__ import annotations

from pathlib import Path

import pytest

import gbdraw.labels.circular_fit as circular_fit_module
from tests.utils.render_stage_equivalence import capture_cli_request, render_captured

TEST_INPUTS = Path(__file__).parent / "test_inputs"


def test_a_refusal_lists_the_same_rows_however_soon_the_bound_comes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # Lambda at 12 pt with gene names first overprints two embedded labels, so
    # the render refuses and searches the font, radial and separate-features rows.
    priority = tmp_path / "priority.tsv"
    priority.write_text("CDS\tgene,product,locus_tag\n", encoding="utf-8")
    captured = capture_cli_request(
        "circular",
        [
            "--gbk", str(TEST_INPUTS / "NC_001416.gb"),
            "--separate_strands", "--track_type", "tuckin", "-l", "left",
            "--labels", "out", "--qualifier_priority", str(priority), "--label_font_size", "12",
            "--center_reserved_radius", "150",
            "-o", str(tmp_path / "out"), "-f", "svg",
        ],
    )

    on_any_host = render_captured(captured, None)
    # A host so slow that the bound has passed before the first row is searched.
    monkeypatch.setattr(circular_fit_module, "SWITCH_SEARCH_SECONDS", 0.0)
    on_a_slow_host = render_captured(captured, None)

    assert on_any_host.error is not None and on_any_host.error[0] == "LabelsDoNotFitError"
    assert "Not searched within" not in on_any_host.error[1]
    assert on_a_slow_host.error == on_any_host.error
