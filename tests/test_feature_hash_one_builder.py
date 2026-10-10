"""One builder of the feature hash (OV-401, OV-412).

A drawn feature's SVG IDs carry the hash of the feature in its source record,
also when its record is cropped or reverse-complemented, so the feature
catalog binds to a plain SVG under the strict validator and every export
writes the same stable feature ID.
"""

from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.api import enrich_svg
from gbdraw.api.options import LinearDiagramOptions, LinearOutputOptions
from gbdraw.api.request_render import render_request
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.exceptions import GbdrawError
from gbdraw.io.regions import parse_region_spec
from gbdraw.linear import linear_main

_SUFFIX_RE = re.compile(r"(__instance_.*|_record_\d+)+$")


def _write_records(tmp_path: Path) -> list[Path]:
    features = [
        SeqFeature(
            SimpleLocation(start, start + 300, strand=strand),
            type="CDS",
            qualifiers={"product": [f"p{start}"], "translation": ["M"]},
        )
        for start, strand in ((100, 1), (900, -1))
    ]
    paths = []
    for name in ("recA", "recB"):
        record = SeqRecord(Seq("ATGC" * 500), id=name, name=name, features=features)
        record.annotations["molecule_type"] = "DNA"
        path = tmp_path / f"{name}.gb"
        SeqIO.write(record, path, "genbank")
        paths.append(path)
    return paths


def _catalog(interactive_svg: str) -> dict:
    root = ET.fromstring(interactive_svg)
    metadata = next(
        element
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "metadata"
        and element.get("id") == "gbdraw-interactive-feature-metadata"
    )
    return json.loads(metadata.text or "{}")


def _render(tmp_path: Path, case: str) -> tuple[str, str]:
    """Return the plain and interactive SVG of recA + a transformed recB."""

    paths = _write_records(tmp_path)
    prefix = tmp_path / "out"
    if case == "cli-reverse":
        linear_main([
            "--gbk", *map(str, paths),
            "--reverse_complement", "0", "--reverse_complement", "1",
            "-f", "svg,interactive_svg", "-o", str(prefix),
        ])
    else:
        region = {
            "crop": "recB:51-1800",
            "crop-reverse": "recB:51-1800:rc",
        }.get(case)
        records = (
            RecordInput(GenBankInputSource(paths[0]), record_key="recA"),
            RecordInput(
                GenBankInputSource(paths[1]),
                record_key="recB",
                region=None if region is None else parse_region_spec(region),
                # An alignment direction choice reverses a record through its
                # presentation (Review alignment options).
                presentation=RecordPresentation(
                    reverse_complement=case == "alignment-reverse"
                ),
            ),
        )
        render_request(
            LinearDiagramRequest(
                records=records,
                options=LinearDiagramOptions(
                    selected_features_set=("CDS",),
                    output=LinearOutputOptions(legend="none"),
                ),
                output=RenderOutputRequest(
                    output_prefix="out",
                    output_directory=tmp_path,
                    formats=("svg", "interactive_svg"),
                ),
            )
        )
    return (
        prefix.with_suffix(".svg").read_text(encoding="utf-8"),
        prefix.with_suffix(".interactive.svg").read_text(encoding="utf-8"),
    )


def _feature_elements(svg: str) -> list[ET.Element]:
    return [
        element
        for element in ET.fromstring(svg).iter()
        if element.get("data-gbdraw-feature-id")
        and element.get("data-gbdraw-auto-feature-underlay") is None
    ]


CASES = ("cli-reverse", "crop", "crop-reverse", "alignment-reverse")


@pytest.mark.parametrize("case", CASES)
def test_transformed_record_draws_source_feature_hash(tmp_path: Path, case: str) -> None:
    plain, interactive = _render(tmp_path, case)
    item = _catalog(interactive)["items"][0]
    record_b = item["recordKeys"][1]
    stable = {
        (row["recordKey"], row["biologicalFeatureId"]): row.get("stableFeatureId")
        or row["biologicalFeatureId"]
        for row in item["biologicalFeatures"]
    }
    rows = [row for row in item["features"] if row["recordKey"] == record_b]
    assert rows

    plain_elements = {element.get("id"): element for element in _feature_elements(plain)}
    for row in rows:
        expected = stable[(row["recordKey"], row["biologicalFeatureId"])]
        # The handle is the source hash plus the record/instance suffixes.
        assert _SUFFIX_RE.sub("", row["svgId"]) == expected
        element = plain_elements[row["svgId"]]
        assert _SUFFIX_RE.sub("", element.get("data-gbdraw-feature-id")) == expected
        # OV-412: the plain SVG (the web app's Result SVG) already writes the
        # stable ID that the interactive export writes.
        assert element.get("data-gbdraw-stable-feature-id") == expected

    interactive_ids = {
        element.get("id"): element.get("data-gbdraw-stable-feature-id")
        for element in _feature_elements(interactive)
    }
    assert interactive_ids == {
        element_id: element.get("data-gbdraw-stable-feature-id")
        for element_id, element in plain_elements.items()
    }

    # OV-401: the strict validator binds gbdraw's own catalog to the plain SVG.
    enrich_svg(
        plain,
        result_index=item["resultIndex"],
        result_name=item["resultName"],
        feature_catalog=_catalog(interactive),
    )


@pytest.mark.parametrize("case", ("plain", "cli-reverse", "crop-reverse"))
def test_catalog_row_swapped_within_a_record_is_rejected(tmp_path: Path, case: str) -> None:
    if case == "plain":
        paths = _write_records(tmp_path)
        prefix = tmp_path / "out"
        linear_main(["--gbk", *map(str, paths), "-f", "svg,interactive_svg", "-o", str(prefix)])
        plain = prefix.with_suffix(".svg").read_text(encoding="utf-8")
        interactive = prefix.with_suffix(".interactive.svg").read_text(encoding="utf-8")
    else:
        plain, interactive = _render(tmp_path, case)
    catalog = _catalog(interactive)
    item = catalog["items"][0]
    rows = [
        row
        for row in item["features"]
        if row["recordKey"] == item["recordKeys"][1]
    ]
    assert len(rows) == 2
    rows[0]["biologicalFeatureId"], rows[1]["biologicalFeatureId"] = (
        rows[1]["biologicalFeatureId"],
        rows[0]["biologicalFeatureId"],
    )
    with pytest.raises(GbdrawError, match="does not agree with rendered SVG ID"):
        enrich_svg(
            plain,
            result_index=item["resultIndex"],
            result_name=item["resultName"],
            feature_catalog=catalog,
        )


@pytest.mark.parametrize("transform", ("reverse", "crop", "crop-reverse"))
def test_source_hash_selector_matches_on_a_transformed_record(
    tmp_path: Path, transform: str
) -> None:
    """`hash=` names a feature by its source-record hash on every record."""
    import pandas as pd

    from gbdraw.features.ids import source_feature_location_parts
    from gbdraw.features.source import build_source_feature_catalog
    from gbdraw.features.visibility import (
        compile_feature_visibility_rules,
        should_render_feature,
    )
    from gbdraw.io.record_select import reverse_records
    from gbdraw.io.regions import apply_region_specs, parse_region_specs

    source = SeqIO.read(_write_records(tmp_path)[1], "genbank")
    catalog = build_source_feature_catalog(source)
    record = source
    if transform != "reverse":
        record = apply_region_specs(
            [record],
            parse_region_specs(["recB:51-1800" + (":rc" if transform == "crop-reverse" else "")]),
        )[0]
    else:
        record = reverse_records([record], True)[0]
    for entry in catalog:
        rules = compile_feature_visibility_rules(pd.DataFrame(
            [["*", "*", "hash", f"^{entry.stable_feature_id}$", "off"]],
            columns=["record_id", "feature_type", "qualifier", "value", "action"],
        ))
        hidden = [
            feature
            for feature in record.features
            if not should_render_feature(
                feature, ["CDS"], feature_visibility_rules=rules, record_id=record.id
            )
        ]
        # The one feature whose source parts are the row's.
        assert [source_feature_location_parts(feature, record) for feature in hidden] == [
            entry.location_parts
        ]
        assert len(hidden) == 1
