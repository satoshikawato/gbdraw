"""A fixed identity colour scale keeps comparison colours and Legend labels across figures.

By default a comparison colours its matches (Linear) or ring hits (Circular)
from the identity threshold to 100%, so the same 95% hit changes colour when
the threshold changes. ``--identity_scale MIN MAX`` (CLI), ``identity_scale``
(Python), and the ``objects.<blast_match|conservation>.identity_scale``
config leaf fix that domain; one owner (``identity_color_domain``) gives it to
the match drawer, the ring drawer, and both Legends.
"""

from __future__ import annotations

import random
import re
import sys
import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.comparisons.identity_scale import (
    CIRCULAR_IDENTITY_SCALE_PATH,
    LINEAR_IDENTITY_SCALE_PATH,
    IdentityColorDomain,
    identity_color_domain,
    identity_scale_from_config,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS

IDENTITIES = (60.0, 80.0, 95.0)


def _record(name: str) -> SeqRecord:
    rng = random.Random(name)
    record = SeqRecord(Seq("".join(rng.choice("ACGT") for _ in range(3000))), id=name, name=name)
    record.annotations["molecule_type"] = "DNA"
    return record


def _hits(query: str, subject: str) -> pd.DataFrame:
    rows = [
        [query, subject, identity, 400, 0, 0, 1 + 800 * index, 400 + 800 * index,
         1 + 800 * index, 400 + 800 * index, 1e-50, 500.0]
        for index, identity in enumerate(IDENTITIES)
    ]
    return pd.DataFrame(rows, columns=COMPARISON_COLUMNS)


def _match_fills(svg: str) -> dict[float, str]:
    fills = {}
    for element in ET.fromstring(svg).iter():
        identity = element.get("data-identity")
        if element.tag.rsplit("}", 1)[-1] == "path" and identity is not None:
            fills[float(identity)] = str(element.get("fill"))
    return fills


def _scale_labels(svg: str) -> list[str]:
    """The distinct percent labels, in drawn order (a Linear SVG may hold a Legend copy)."""
    return list(dict.fromkeys(re.findall(r">(\d+(?:\.\d+)?%)<", svg)))


def _run_cli(tmp_path, monkeypatch, mode: str, *args: str) -> str:
    from gbdraw import cli

    output = tmp_path / f"out-{abs(hash(args))}"
    monkeypatch.setattr(sys, "argv", ["gbdraw", mode, *args, "-o", str(output), "-f", "svg"])
    cli.main()
    return output.with_suffix(".svg").read_text()


# The owner


def test_the_default_domain_runs_from_the_threshold_to_100() -> None:
    assert identity_color_domain(70, None) == IdentityColorDomain(70.0, 100.0)
    assert identity_color_domain(70, (40.0, 90.0)) == IdentityColorDomain(40.0, 90.0)


def test_the_factor_is_clamped_and_a_domain_without_width_colours_at_its_maximum() -> None:
    domain = IdentityColorDomain(40.0, 90.0)
    assert domain.factor(65) == pytest.approx(0.5)
    assert domain.factor(10) == 0.0
    assert domain.factor(99) == 1.0
    assert IdentityColorDomain(100.0, 100.0).factor(100) == 1.0
    assert IdentityColorDomain(100.0, 100.0).factor(99.9) == 0.0


@pytest.mark.parametrize("value", [None, "auto", "AUTO", ""])
def test_auto_is_no_fixed_scale(value: object) -> None:
    assert identity_scale_from_config(value, path=LINEAR_IDENTITY_SCALE_PATH) is None


@pytest.mark.parametrize("value", [[50, 100], (0, 12.5)])
def test_two_ordered_percentages_are_a_fixed_scale(value: object) -> None:
    assert identity_scale_from_config(value, path=LINEAR_IDENTITY_SCALE_PATH) == tuple(float(item) for item in value)  # type: ignore[union-attr]


@pytest.mark.parametrize(
    ("value", "reason"),
    [
        ([90, 50], "INCREASING"),
        ([50, 50], "INCREASING"),
        ([100, 100], "INCREASING"),
        ([-1, 50], "PERCENT"),
        ([50, 101], "PERCENT"),
        ([50], "FINITE"),
        ("50 100", "FINITE"),
        ([True, 100], "FINITE"),
        ([50, float("nan")], "FINITE"),
    ],
)
def test_other_values_are_refused_with_the_config_path_and_reason(value: object, reason: str) -> None:
    with pytest.raises(ValidationError) as error:
        apply_config_overrides(None, {CIRCULAR_IDENTITY_SCALE_PATH: value})
    assert error.value.diagnostic == {"code": "INPUT_INVALID", "reason": reason, "configPath": CIRCULAR_IDENTITY_SCALE_PATH}


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("ends", [("90", "40"), ("100", "100"), ("50", "101")])
def test_the_cli_refuses_a_scale_naming_the_option(mode: str, ends: tuple[str, str], capsys) -> None:
    from gbdraw import circular, linear

    get_args = circular._get_args if mode == "circular" else linear._get_args
    with pytest.raises(SystemExit) as raised:
        get_args(["--gbk", "dummy.gb", "--identity_scale", *ends])
    assert raised.value.code == 2
    assert "argument --identity_scale: must be two percent identities MIN MAX" in capsys.readouterr().err


def test_the_packaged_config_has_no_fixed_scale() -> None:
    config = apply_config_overrides(None, None)
    assert config.objects.blast_match.identity_scale is None
    assert config.objects.conservation.identity_scale is None
    assert apply_config_overrides(None, {LINEAR_IDENTITY_SCALE_PATH: [50, 100]}).objects.blast_match.identity_scale == (
        50.0, 100.0,
    )


# Linear: the match drawer and the Legend


def _linear_inputs(tmp_path) -> list[str]:
    paths = []
    for name in ("q1", "s1"):
        path = tmp_path / f"{name}.gb"
        SeqIO.write([_record(name)], str(path), "genbank")
        paths.append(str(path))
    blast = tmp_path / "q1_s1.tsv"
    _hits("q1", "s1").to_csv(blast, sep="\t", header=False, index=False)
    return [*paths, str(blast)]


def test_linear_cli_fixed_scale_keeps_colours_and_legend_when_the_threshold_changes(tmp_path, monkeypatch) -> None:
    query, subject, blast = _linear_inputs(tmp_path)
    runs = {
        threshold: _run_cli(
            tmp_path, monkeypatch, "linear", "--gbk", query, subject, "-b", blast,
            "--identity", threshold, "--identity_scale", "40", "100",
        )
        for threshold in ("50", "70")
    }

    low, high = _match_fills(runs["50"]), _match_fills(runs["70"])
    assert set(low) == {60.0, 80.0, 95.0} and set(high) == {80.0, 95.0}
    assert low[80.0] == high[80.0] and low[95.0] == high[95.0]
    assert _scale_labels(runs["50"]) == _scale_labels(runs["70"]) == ["40%", "100%"]
    factors = [float(value) for value in re.findall(r'data-identity-factor="([^"]+)"', runs["70"])]
    assert factors == pytest.approx([(80 - 40) / 60, (95 - 40) / 60])


def test_linear_cli_default_scale_still_follows_the_threshold(tmp_path, monkeypatch) -> None:
    query, subject, blast = _linear_inputs(tmp_path)
    low = _run_cli(tmp_path, monkeypatch, "linear", "--gbk", query, subject, "-b", blast, "--identity", "50")
    high = _run_cli(tmp_path, monkeypatch, "linear", "--gbk", query, subject, "-b", blast, "--identity", "70")

    assert _scale_labels(low) == ["50%", "100%"] and _scale_labels(high) == ["70%", "100%"]
    assert _match_fills(low)[95.0] != _match_fills(high)[95.0]


def test_linear_python_option_matches_the_cli(tmp_path, monkeypatch) -> None:
    import gbdraw

    query, subject, blast = _linear_inputs(tmp_path)
    cli_svg = _run_cli(
        tmp_path, monkeypatch, "linear", "--gbk", query, subject, "-b", blast,
        "--identity", "50", "--identity_scale", "40", "90",
    )
    diagram = gbdraw.draw_linear(
        [_record("q1"), _record("s1")],
        options=gbdraw.LinearOptions(
            comparisons=gbdraw.LinearComparisonOptions(blast_files=[blast], identity_scale=(40, 90)),
            thresholds=gbdraw.Thresholds(identity=50),
        ),
    )

    assert _match_fills(diagram.to_svg()) == _match_fills(cli_svg)
    assert _scale_labels(diagram.to_svg()) == _scale_labels(cli_svg) == ["40%", "90%"]


def test_python_option_refuses_a_reversed_scale() -> None:
    import gbdraw

    with pytest.raises(ValidationError) as linear_error:
        gbdraw.LinearComparisonOptions(identity_scale=(90, 40))
    with pytest.raises(ValidationError) as circular_error:
        gbdraw.ComparisonRingOptions(identity_scale=(90, 40))
    assert linear_error.value.diagnostic["configPath"] == LINEAR_IDENTITY_SCALE_PATH
    assert circular_error.value.diagnostic["configPath"] == CIRCULAR_IDENTITY_SCALE_PATH


# Circular: the ring drawer and the Legend


def _ring_fills(svg: str) -> dict[float, str]:
    fills = {}
    for element in ET.fromstring(svg).iter():
        identity = element.get("data-identity")
        if identity is not None and element.get("fill"):
            fills[float(identity)] = str(element.get("fill"))
    return fills


def _circular_inputs(tmp_path) -> list[str]:
    gbk = tmp_path / "ref.gb"
    SeqIO.write([_record("ref")], str(gbk), "genbank")
    blast = tmp_path / "ring.tsv"
    _hits("other", "ref").to_csv(blast, sep="\t", header=False, index=False)
    return [str(gbk), str(blast)]


def test_circular_cli_fixed_scale_keeps_ring_colours_and_legend_when_the_threshold_changes(
    tmp_path, monkeypatch,
) -> None:
    gbk, blast = _circular_inputs(tmp_path)
    runs = {
        threshold: _run_cli(
            tmp_path, monkeypatch, "circular", "--gbk", gbk, "--conservation_blast", blast,
            "--conservation_reference", "subject", "--bitscore", "0",
            "--identity", threshold, "--identity_scale", "40", "100",
        )
        for threshold in ("50", "70")
    }

    low, high = _ring_fills(runs["50"]), _ring_fills(runs["70"])
    assert set(low) == {60.0, 80.0, 95.0} and set(high) == {80.0, 95.0}
    assert low[80.0] == high[80.0] and low[95.0] == high[95.0]
    assert _scale_labels(runs["50"]) == _scale_labels(runs["70"]) == ["40%", "100%"]


def test_circular_python_option_matches_the_cli(tmp_path, monkeypatch) -> None:
    import gbdraw

    gbk, blast = _circular_inputs(tmp_path)
    cli_svg = _run_cli(
        tmp_path, monkeypatch, "circular", "--gbk", gbk, "--conservation_blast", blast,
        "--conservation_reference", "subject", "--bitscore", "0", "--identity", "50", "--identity_scale", "30", "90",
    )
    rings = gbdraw.ComparisonRingOptions(
        tracks=[gbdraw.ComparisonRingTrackOptions(source=blast)], reference="subject", identity_scale=(30, 90),
    )
    diagram = gbdraw.draw_circular(
        _record("ref"),
        options=gbdraw.CircularOptions(comparison_rings=rings, thresholds=gbdraw.Thresholds(identity=50, bitscore=0)),
    )

    assert _ring_fills(diagram.to_svg()) == _ring_fills(cli_svg)
    assert _scale_labels(diagram.to_svg()) == _scale_labels(cli_svg) == ["30%", "90%"]
