"""Compare an actual browser download with its effective Web draft rows."""

from dataclasses import asdict
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from gbdraw.annotations import CoordinateSpan, read_annotation_table


def camel_fields(value):
    if isinstance(value, dict):
        return {
            key.split("_")[0] + "".join(part.title() for part in key.split("_")[1:]): camel_fields(item)
            for key, item in value.items()
        }
    if isinstance(value, (tuple, list)):
        return [camel_fields(item) for item in value]
    return value


rows = []
for annotation_set in read_annotation_table(sys.argv[1]):
    for item in annotation_set.annotations:
        target = camel_fields(asdict(item.target))
        target["kind"] = "coordinateSpan" if isinstance(item.target, CoordinateSpan) else "featureSpan"
        record = item.target.record
        target["record"] = (
            None if record is None else
            {"kind": "recordIndex", "index": record.record_index} if record.record_index is not None else
            {"kind": "recordId", "value": record.record_id}
        )
        rows.append({
            "setId": annotation_set.id, "id": item.id, "target": target,
            "mark": item.mark, "label": item.label or "", "lane": item.lane,
            "legendLabel": item.legend_label or annotation_set.legend_label,
            "style": camel_fields(asdict(item.style or annotation_set.default_style)),
        })
expected = json.loads(Path(sys.argv[2]).read_text(encoding="utf-8"))
assert rows == expected, json.dumps({"python": rows, "web": expected}, indent=2)
print(f"Python reader: {len(rows)} effective annotation rows match the downloaded Web draft")
