"""Read-only probes of issue #600; run from the repository root with PYTHONPATH=."""
from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations import (
    AnnotationSet, FeatureSelector, FeatureSpan, RegionAnnotation,
    annotation_sets_from_dataframe, resolve_annotations,
)
from gbdraw.exceptions import ValidationError
from gbdraw.legend.table import prepare_legend_table
from gbdraw.tracks.circular import parse_circular_track_slot


def observe(name, operation):
    try:
        value = operation()
        return {"case": name, "value": value}
    except (ValidationError, ValueError) as error:
        return {"case": name, "error_type": type(error).__name__, "error": str(error)}


record = SeqRecord(Seq("A" * 100), id="r1")
record.features = [SeqFeature(FeatureLocation(10, 20), type="CDS", qualifiers={"gene": ["present"]})]


def annotation_probe(*genes):
    target = FeatureSpan(None, tuple(FeatureSelector(gene, "gene") for gene in genes))
    bundle = resolve_annotations((AnnotationSet("s", (RegionAnnotation("a", target),)),), [record], mode="linear")
    return {"count": len(bundle.annotations), "warnings": [item.code for item in bundle.warnings]}


def legend_probe():
    rules = pd.DataFrame([
        ["CDS", "gene", "abc", "#112233", "Transporter"],
        ["CDS", "gene", "mfs", "#445566", "Transporter"],
    ], columns=["feature_type", "qualifier_key", "value", "color", "caption"])
    gc = SimpleNamespace(stroke_color="#000000", stroke_width=1, high_fill_color="#aaaaaa", low_fill_color="#bbbbbb", dinucleotide="GC")
    skew = SimpleNamespace(stroke_color="#000000", stroke_width=1, high_fill_color="#cccccc", low_fill_color="#dddddd")
    features = SimpleNamespace(color_table=rules, default_colors=pd.DataFrame([["CDS", "#cccccc"]], columns=["feature_type", "color"]), block_stroke_color="#000000", block_stroke_width=1)
    result = prepare_legend_table(gc, skew, features, ["CDS"], used_color_rules={("Transporter", "#112233"), ("Transporter", "#445566")}, show_gc=False, show_skew=False, show_depth=False)
    return {caption: entry["fill"] for caption, entry in result.items()}


observations = [
    observe("BUG-08 auxiliary column", lambda: len(annotation_sets_from_dataframe(pd.DataFrame([{"set_id": "s", "id": "a", "mark": "band", "start": 1, "end": 5, "notes": "auxiliary"}])))),
    observe("BUG-09 all missing", lambda: annotation_probe("absent")),
    observe("BUG-09 partial match", lambda: annotation_probe("present", "absent")),
    observe("BUG-09 control", lambda: annotation_probe("present")),
    observe("BUG-10 two used colors", legend_probe),
    observe("BUG-14 circular CLI slot", lambda: parse_circular_track_slot("gc:dinucleotide_content@inner_gap_px=10px").inner_gap_px),
]
print(json.dumps(observations, indent=2, ensure_ascii=False))
