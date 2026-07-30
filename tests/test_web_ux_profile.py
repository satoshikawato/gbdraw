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
