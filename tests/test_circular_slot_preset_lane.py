"""OV-1326: an explicit Circular features slot's ``side`` sets its lane; the track
preset's ``lane_direction`` never conflicts with it. Two explicit values that
conflict are still refused where they are entered."""

from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw import CircularOptions, CircularTrackOptions, draw_circular, read_genbank
from gbdraw.exceptions import ValidationError
from gbdraw.interface import AnnotationOptions

HUMAN = Path(__file__).resolve().parents[1] / "gbdraw" / "web" / "tutorial-data" / "human-mitochondrion" / "HmmtDNA.gbk"
ANNOTATIONS = (
    "set_id\tid\tmark\tstart\tend\tlabel\tfill\tfill_opacity\tlegend_label\n"
    "review\twindow\tband\t1000\t5000\tReview window\t#f59e0b\t0.25\tReview region\n"
)
# docs/RECIPES.md "Highlight named regions"; the Tutorials add `,lane_direction=split`.
RECIPE_SLOTS = (
    "review:annotations@set_id=review,side=outside,w=28px",
    "features:features@side=overlay",
    "ticks:ticks@side=inside",
)
LANE_FOR_SIDE = {"inside": "inside", "outside": "outside", "overlay": "split"}


def _with_lane(slots: tuple[str, ...], lane: str) -> tuple[str, ...]:
    return tuple(f"{slot},lane_direction={lane}" if slot.startswith("features:") else slot for slot in slots)


@pytest.mark.circular
def test_the_recipe_command_draws_what_the_tutorial_spelling_draws(tmp_path: Path, gbdraw_runner) -> None:
    table = tmp_path / "annotations.tsv"
    table.write_text(ANNOTATIONS, encoding="utf-8")
    svgs = []
    for name, slots in (("recipe", RECIPE_SLOTS), ("tutorial", _with_lane(RECIPE_SLOTS, "split"))):
        extra = ["--annotation_table", str(table)]
        for slot in slots:
            extra += ["--circular_track_slot", slot]
        returncode, output, svg_path = gbdraw_runner.run("circular", [HUMAN], name, tmp_path, extra_args=extra)
        assert returncode == 0, output
        assert "Traceback" not in output
        svgs.append(svg_path.read_text(encoding="utf-8").replace(name, "OUT"))
    assert svgs[0] == svgs[1]


def test_the_python_recipe_draws_what_the_tutorial_spelling_draws(tmp_path: Path) -> None:
    table = tmp_path / "annotations.tsv"
    table.write_text(ANNOTATIONS, encoding="utf-8")
    records = read_genbank(str(HUMAN))
    svgs = [
        draw_circular(
            records,
            options=CircularOptions(
                annotations=AnnotationOptions(table_file=str(table)),
                tracks=CircularTrackOptions(slots=slots),
            ),
        ).to_svg()
        for slots in (RECIPE_SLOTS, _with_lane(RECIPE_SLOTS, "split"))
    ]
    assert svgs[0] == svgs[1]


@pytest.mark.parametrize("track_type", ["tuckin", "middle", "spreadout"])
@pytest.mark.parametrize("side", ["inside", "outside", "overlay"])
def test_an_explicit_side_draws_the_lane_it_implies_under_every_preset(track_type: str, side: str) -> None:
    records = read_genbank(str(HUMAN))
    svgs = [
        draw_circular(
            records,
            options=CircularOptions(
                config_overrides={"canvas.circular.track_type": track_type},
                tracks=CircularTrackOptions(slots=slots),
            ),
        ).to_svg()
        for slots in ((f"features:features@side={side}",), (f"features:features@side={side},lane_direction={LANE_FOR_SIDE[side]}",))
    ]
    assert svgs[0] == svgs[1]


def test_two_explicit_values_that_conflict_are_still_refused(capsys: pytest.CaptureFixture[str]) -> None:
    import gbdraw.circular as circular_cli_module

    slot = "features:features@side=overlay,lane_direction=inside"
    with pytest.raises(ValidationError, match=r"features slot 'features' has conflicting side='overlay' and lane_direction='inside'"):
        CircularTrackOptions(slots=(slot,))
    with pytest.raises(SystemExit) as excinfo:
        circular_cli_module._get_args(["--gbk", str(HUMAN), "--circular_track_slot", slot])
    assert excinfo.value.code == 2
    assert "conflicting side='overlay' and lane_direction='inside'" in capsys.readouterr().err
