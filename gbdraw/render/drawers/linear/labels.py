#!/usr/bin/env python
# coding: utf-8

from ....svg.text_path import generate_text_path


class LabelDrawer:
    def __init__(self, config_dict: dict) -> None:
        self.config_dict = config_dict

    def add_label(self, label_entry, group):
        feature_label_text = label_entry["label_text"]
        middle_x = label_entry["middle_x"]
        middle_y = label_entry["middle_y"]
        rotation_deg = float(label_entry.get("rotation_deg", 0.0))
        label_path = generate_text_path(
            text=feature_label_text,
            title_x=middle_x,
            title_y=middle_y,
            interval=0,
            font_size=label_entry["font_size"],
            font_weight="normal",
            font=label_entry["font_family"],
            dominant_baseline="central",
            text_anchor=label_entry.get("text_anchor", "middle"),
        )
        if rotation_deg != 0.0:
            label_path.rotate(rotation_deg, center=(middle_x, middle_y))
        group.add(label_path)
        return group

    def draw(self, label_entry, group):
        return self.add_label(label_entry, group)


__all__ = ["LabelDrawer"]


