"""The mode-profile data the Web app reads, generated into ``mode-profiles.generated.js``.

``tools/generate_mode_profiles.py`` writes it; each value comes from its owner module.
"""

from __future__ import annotations

from typing import get_args

from gbdraw.api.region_drawing import NotCarriedKind, StartAfterEndReason
from gbdraw.auto_sizes import (
    ADAPTABLE_SETTINGS,
    linear_label_size_records,
    scalar_axis_tick_font_rule,
    tick_and_window_tiers,
)
from gbdraw.comparisons.identity_scale import IDENTITY_SCALE_PATHS
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.labels import LinearLabelScope
from gbdraw.config.toml import load_config_toml
from gbdraw.core.figure_length import figure_length_unit_table
from gbdraw.diagrams.circular.presets import circular_auto_slot_sizes, circular_track_preset_table
from gbdraw.labels.label_style import (
    DEFAULT_ITALIC_GENE_LABELS,
    DEFAULT_LABEL_LINE,
    DEFAULT_LABEL_PROTECTION,
    LABEL_LINES,
    LABEL_PROTECTIONS,
)
from gbdraw.mode_profiles import (
    COMPARISON_THRESHOLD_DOMAINS,
    DEFAULT_FEATURE_TYPES,
    MODE_PROFILE_VERSION,
    MODE_PROFILES,
)


def _camel_threshold_name(name: str) -> str:
    head, *tail = name.split("_")
    return head + "".join(part.title() for part in tail)


def mode_profiles_payload() -> dict[str, object]:
    """Return the JSON-compatible profile representation used by the web app."""

    packaged = load_config_toml("gbdraw.data", "config.toml")
    packaged_cfg = GbdrawConfig.from_dict(packaged)
    length_threshold = packaged["labels"]["length_threshold"]
    return {
        "version": MODE_PROFILE_VERSION,
        # The packaged text font list: a resolved config that holds it names no font.
        "defaultFontFamily": str(packaged["objects"]["text"]["font_family"]),
        # The packaged short/long boundary of each mode, the records that set
        # the Auto Linear label font for each label display scope (rules:
        # gbdraw/auto_sizes.py), and the Auto sizes of a Circular row
        # (gbdraw/diagrams/circular/presets.py), the settings of each mode
        # that can return to Auto, by setting ID, and the Auto tick font of a
        # GC content or Depth axis by track thickness.
        "autoSizes": {
            "adaptableSettings": {mode: list(settings) for mode, settings in ADAPTABLE_SETTINGS.items()},
            "circularSlots": circular_auto_slot_sizes(packaged_cfg),
            "linearLabelSizeRecords": {
                scope: linear_label_size_records(scope) for scope in get_args(LinearLabelScope)
            },
            "scalarAxisTickFont": scalar_axis_tick_font_rule(),
            "sizeClassThresholdBp": {mode: int(length_threshold[mode]) for mode in MODE_PROFILES},
            **tick_and_window_tiers(packaged_cfg),
        },
        # The Circular track presets in menu order and their feature lane directions.
        "circularTrackPresets": circular_track_preset_table(),
        "comparisonDomains": {
            _camel_threshold_name(name): {
                "minimum": domain.minimum,
                "maximum": domain.maximum,
                "integer": domain.integer,
                "reason": domain.reason,
            }
            for name, domain in COMPARISON_THRESHOLD_DOMAINS.items()
        },
        "featureTypes": list(DEFAULT_FEATURE_TYPES),
        # The units of canvas.figure_width and canvas.figure_height (gbdraw/core/figure_length.py): the Figure size
        # field's unit choices and the size line's conversion and decimals.
        "figureLength": figure_length_unit_table(),
        # The vocabularies of a region drawing's reply (gbdraw/api/region_drawing.py):
        # the kinds of a part the new drawing does not carry and the reasons a
        # region is refused.
        "regionDrawing": {
            "notCarriedKinds": list(get_args(NotCarriedKind)),
            "refusalReasons": list(get_args(StartAfterEndReason)),
        },
        # The values of labels.<mode>.line and .protection.
        "labelStyleChoices": {"line": list(LABEL_LINES), "protection": list(LABEL_PROTECTIONS)},
        "modes": {
            mode: {
                "comparison": {
                    "evalue": profile.comparison.evalue,
                    "bitscore": profile.comparison.bitscore,
                    "identity": profile.comparison.identity,
                    "alignmentLength": profile.comparison.alignment_length,
                },
                "tracks": {
                    "gc": profile.show_gc,
                    "skew": profile.show_skew,
                },
                "linearAxisColor": profile.linear_axis_color,
                "linearRulerAxisColor": profile.linear_ruler_axis_color,
                # The packaged feature label style of the mode (labels.<mode>.line, .protection,
                # .italic_gene_labels).
                "labelStyle": {
                    "line": DEFAULT_LABEL_LINE[mode],
                    "protection": DEFAULT_LABEL_PROTECTION[mode],
                    "italicGeneLabels": DEFAULT_ITALIC_GENE_LABELS[mode],
                },
                # The config leaf of the mode's fixed identity colour scale (an ``[MIN, MAX]`` override).
                "identityScalePath": IDENTITY_SCALE_PATHS[mode],
            }
            for mode, profile in MODE_PROFILES.items()
        },
    }


__all__ = ["mode_profiles_payload"]
